module FreeRun

using MarbleRun
using MarbleRun: TOML, TOML_Phase
using Marble
using GeometricObjects

using Serialization

function preprocess_input!(input::TOML)
end

function initialize(input::TOML)
    GridState = @NamedTuple begin
        m           :: Float64
        m′          :: Float64
        v           :: Vec{2, Float64}
        v_n         :: Vec{2, Float64}
        m_contacted :: Float64
        vᵣ          :: Vec{2, Float64}
        fc          :: Vec{2, Float64}
        d           :: Vec{2, Float64}
        μ           :: Vec{2, Float64} # [μ, c]
        poly_coef   :: Vec{3, Float64}
        poly_mat    :: Mat{3, 3, Float64, 9}
    end
    L = isa(input.General.interpolation, LinearWLS) ? 3 : 2
    PointState = @NamedTuple begin
        m        :: Float64
        V        :: Float64
        x        :: Vec{2, Float64}
        x0       :: Vec{2, Float64}
        v        :: Vec{2, Float64}
        b        :: Vec{2, Float64}
        fc       :: Vec{2, Float64}
        σ        :: SymmetricSecondOrderTensor{3, Float64, 6}
        ϵ        :: SymmetricSecondOrderTensor{3, Float64, 6}
        ∇v       :: SecondOrderTensor{3, Float64, 9}
        P        :: Float64
        C        :: Mat{2, L, Float64, 2*L}
        r        :: Vec{2, Float64}
        index    :: Int
        matindex :: Int
    end

    coordinate_system = input.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = input.General.domain
    dx = input.General.grid_space
    g = input.General.gravity

    materials = input.Material
    rigidbodies = map(x -> x.model, input.RigidBody)

    grid = Grid(xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    gridstate = generate_gridstate(GridState, grid)
    pointstate = generate_pointstate(PointState, grid, input) do pointstate, matindex
        MarbleRun.initialize_pointstate!(pointstate, materials[matindex], g)
        @. pointstate.matindex = matindex
    end
    t = 0.0

    for dirichlet in input.BoundaryCondition.Dirichlet
        node_indices = findall(x -> dirichlet.region(x[1], x[2]), grid)
        dirichlet.node_indices = node_indices
    end

    t, grid, gridstate, pointstate, rigidbodies
end

function main(input::TOML, phase::TOML_Phase, t, grid::Grid, gridstate::AbstractArray, pointstate::AbstractVector, rigidbodies)

    # General/Output
    dx = input.General.grid_space
    t_start = t
    t_stop = phase.time_stop
    t_step = input.Output.time_interval

    # Material models
    matmodels = map(x -> x.model, input.Material)

    ################
    # Output files #
    ################

    outdir = input.Output.directory
    outputs = Dict{String, Any}()
    if input.Paraview.output
        mkpath(joinpath(outdir, "paraview"))
        outputs["paraview_file"] = joinpath(outdir, "paraview", "output")
        paraview_collection(vtk_save, outputs["paraview_file"])
    end
    if input.Output.snapshots || input.Output.snapshot_last
        mkpath(joinpath(outdir, "snapshots"))
    end
    if any(d -> d.output, input.BoundaryCondition.Dirichlet)
        Dirichlet = input.BoundaryCondition.Dirichlet
        for i in eachindex(Dirichlet)
            if Dirichlet[i].output
                dir = joinpath(outdir, "dirichlet", "$i")
                mkpath(dir)
                history_file = joinpath(dir, "history.csv")
                write(history_file, "disp,force\n")
            end
        end
    end
    if any(d -> d.output, input.RigidBody)
        RigidBody = input.RigidBody
        for i in eachindex(RigidBody)
            if RigidBody[i].output
                dim = ndims(grid)
                dir = joinpath(outdir, "rigidbodies", "$i")
                mkpath(dir)
                history_file = joinpath(dir, "history.csv")
                header = [
                    "t"
                    ["x_$i" for i in 1:dim]
                    ["v_$i" for i in 1:dim]
                    ["ω_$i" for i in 1:3]
                    ["attitude_$i" for i in 1:3]
                ]
                write(history_file, join(header, ",") * "\n")
            end
        end
    end
    if isdefined(input.Injection, :main_output)
        Base.invokelatest(input.Injection.main_output_initialize, (;
            input,
            t,
            grid,
            gridstate,
            pointstate,
            rigidbodies,
        ))
    end

    ##################
    # Run simulation #
    ##################

    space = MPSpace(input.General.interpolation, grid, pointstate.x)
    logger = Logger(t_start, t_stop, t_step; input.General.showprogress)
    update!(logger, t)
    writeoutput(outputs, input, grid, gridstate, pointstate, rigidbodies, t, logindex(logger))

    try
        while !isfinised(logger, t)
            dt = phase.CFL * MarbleRun.safe_minimum(pointstate) do pt
                MarbleRun.timestep(matmodels[pt.matindex], pt, dx)
            end
            MarbleRun.advancestep!(grid, gridstate, pointstate, rigidbodies, space, dt, input, phase)

            if input.Output.quickview
                update!(logger, t += dt; print = MarbleRun.quickview_sparsity_pattern(space.spat))
            else
                update!(logger, t += dt)
            end

            if islogpoint(logger)
                if input.Advanced.reorder_pointstate
                    Marble.reorder_pointstate!(pointstate, space)
                end
                writeoutput(outputs, input, grid, gridstate, pointstate, rigidbodies, t, logindex(logger))
            end
        end
    catch e
        writeoutput(outputs, input, grid, gridstate, pointstate, rigidbodies, t, "error")
        rethrow()
    end

    if input.Output.snapshot_last
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot_last"),
            (; t, grid, gridstate, pointstate, rigidbodies)
        )
    end

    t
end

function writeoutput(
        outputs::Dict{String, Any},
        input::TOML,
        grid::Grid,
        gridstate::AbstractArray,
        pointstate::AbstractVector,
        rigidbodies::Vector,
        t::Real,
        output_index,
    )
    if input.Paraview.output
        compress = true
        paraview_file = outputs["paraview_file"]
        paraview_collection(paraview_file, append = true) do pvd
            vtk_multiblock(string(paraview_file, output_index)) do vtm
                vtk_grid(vtm, pointstate.x; compress) do vtk
                    MarbleRun.writevtk(vtk, input.Paraview.PointState, pointstate)
                end
                for rigidbody in rigidbodies
                    vtk_grid(vtm, rigidbody)
                end
                if input.Paraview.GridState !== nothing
                    vtk_grid(vtm, grid; compress) do vtk
                        MarbleRun.writevtk(vtk, input.Paraview.GridState, gridstate)
                    end
                end
                pvd[t] = vtm
            end
        end
    end

    if input.Output.snapshots
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot$output_index"),
            (; t, grid, gridstate, pointstate, rigidbodies)
        )
    end

    if any(d -> d.output, input.BoundaryCondition.Dirichlet)
        Dirichlet = input.BoundaryCondition.Dirichlet
        for i in eachindex(Dirichlet)
            dirichlet = Dirichlet[i]
            if dirichlet.output
                history_file = joinpath(input.Output.directory, "dirichlet", "$i", "history.csv")
                open(history_file, "a") do io
                    disp = dirichlet.displacement
                    force = dirichlet.reaction_force
                    write(io, join([disp, force], ",") * "\n")
                end
            end
        end
    end

    if any(d -> d.output, input.RigidBody)
        RigidBody = input.RigidBody
        for i in eachindex(RigidBody)
            if RigidBody[i].output
                history_file = joinpath(input.Output.directory, "rigidbodies", "$i", "history.csv")
                open(history_file, "a") do io
                    values = [
                        t
                        centroid(rigidbodies[i])
                        rigidbodies[i].v
                        rigidbodies[i].ω
                        attitude(rigidbodies[i])
                    ]
                    write(io, join(values, ",") * "\n")
                end
            end
        end
    end

    if isdefined(input.Injection, :main_output)
        Base.invokelatest(input.Injection.main_output, (;
            input,
            grid,
            gridstate,
            pointstate,
            rigidbodies,
            t,
            output_index,
        ))
    end
end

end

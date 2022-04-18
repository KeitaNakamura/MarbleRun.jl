module FreeRun

using PoingrSimulator
using PoingrSimulator: Input, Input_Phase
using Poingr
using GeometricObjects

using Serialization

function preprocess_input!(input::Input)
    for mat in input.Material
        if hasproperty(mat, :friction_with_rigidbodies)
            @assert length(mat.friction_with_rigidbodies) == length(input.RigidBody)
        end
    end
end

function initialize(input::Input)
    NodeState = @NamedTuple begin
        m::Float64
        m′::Float64
        v::Vec{2, Float64}
        v_n::Vec{2, Float64}
        m_contacted::Float64
        vᵣ::Vec{2, Float64}
        fc::Vec{2, Float64}
        d::Vec{2, Float64}
        μ::Vec{2, Float64} # [μ, c]
    end
    L = isa(input.General.transfer, LinearWLS) ? 3 : 2
    PointState = @NamedTuple begin
        m::Float64
        V::Float64
        x::Vec{2, Float64}
        x0::Vec{2, Float64}
        v::Vec{2, Float64}
        b::Vec{2, Float64}
        fc::Vec{2, Float64}
        σ::SymmetricSecondOrderTensor{3, Float64, 6}
        ϵ::SymmetricSecondOrderTensor{3, Float64, 6}
        ∇v::SecondOrderTensor{3, Float64, 9}
        C::Mat{2, L, Float64, 2*L}
        r::Vec{2, Float64}
        index::Int
        matindex::Int
    end

    coordinate_system = input.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = input.General.domain
    dx = input.General.grid_space
    g = input.General.gravity

    materials = input.Material
    rigidbodies = map(x -> x.model, input.RigidBody)

    grid = Grid(NodeState, input.General.interpolation, xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate(PointState, grid, input) do pointstate, matindex
        PoingrSimulator.initialize_pointstate!(pointstate, materials[matindex], g)
        @. pointstate.matindex = matindex
    end
    t = 0.0

    for dirichlet in input.BoundaryCondition.Dirichlet
        active_nodes = map(x -> dirichlet.inbounds(x[1], x[2]), grid)
        setbounds!(grid, active_nodes)
        dirichlet.active_nodes = active_nodes
    end

    t, grid, pointstate, rigidbodies
end

function main(input::Input, phase::Input_Phase, t, grid::Grid{dim}, pointstate, rigidbodies) where {dim}

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
    if input.Output.paraview
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
        input.Injection.main_output_initialize((;
            input,
            t,
            grid,
            pointstate,
            rigidbodies,
        ))
    end

    ##################
    # Run simulation #
    ##################

    cache = MPCache(grid, pointstate.x)
    logger = Logger(t_start, t_stop, t_step; input.General.showprogress)
    update!(logger, t)
    writeoutput(outputs, input, grid, pointstate, rigidbodies, t, logindex(logger))

    try
        while !isfinised(logger, t)
            dt = phase.CFL * PoingrSimulator.safe_minimum(pointstate) do pt
                PoingrSimulator.timestep(matmodels[pt.matindex], pt, dx)
            end
            PoingrSimulator.advancestep!(grid, pointstate, rigidbodies, cache, dt, input, phase)

            if input.Output.quickview
                update!(logger, t += dt; print = PoingrSimulator.quickview_sparsity_pattern(cache.spat))
            else
                update!(logger, t += dt)
            end

            if islogpoint(logger)
                if input.Advanced.reorder_pointstate
                    Poingr.reorder_pointstate!(pointstate, cache)
                end
                writeoutput(outputs, input, grid, pointstate, rigidbodies, t, logindex(logger))
            end
        end
    catch e
        writeoutput(outputs, input, grid, pointstate, rigidbodies, t, "error")
        rethrow()
    end

    if input.Output.snapshot_last
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot_last"),
            (; t, grid, pointstate, rigidbodies)
        )
    end

    t
end

function writeoutput(
        outputs::Dict{String, Any},
        input::Input,
        grid::Grid,
        pointstate::AbstractVector,
        rigidbodies::Vector,
        t::Real,
        output_index,
    )
    if input.Output.paraview
        compress = true
        paraview_file = outputs["paraview_file"]
        paraview_collection(paraview_file, append = true) do pvd
            vtk_multiblock(string(paraview_file, output_index)) do vtm
                vtk_points(vtm, pointstate.x; compress) do vtk
                    PoingrSimulator.write_vtk_points(vtk, pointstate)
                end
                for rigidbody in rigidbodies
                    vtk_grid(vtm, rigidbody)
                end
                if input.Output.paraview_grid
                    vtk_grid(vtm, grid; compress) do vtk
                        vtk["nodal contact force"] = vec(grid.state.fc)
                        vtk["nodal contact distance"] = vec(grid.state.d)
                        vtk["nodal friction"] = vec(grid.state.μ[1])
                    end
                end
                pvd[t] = vtm
            end
        end
    end

    if input.Output.snapshots
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot$output_index"),
            (; t, grid, pointstate, rigidbodies)
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
        input.Injection.main_output((;
            input,
            grid,
            pointstate,
            rigidbodies,
            t,
            output_index,
        ))
    end
end

end

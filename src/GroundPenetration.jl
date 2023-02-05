module GroundPenetration

using MarbleRun
using MarbleRun: Input, Input_Phase
using GeometricObjects

using Serialization

function preprocess_input!(input::Input)
    input.Material = input.SoilLayer
    input.BoundaryCondition.sides = [
        "-x" => CoulombFriction(:slip),
        "+x" => CoulombFriction(:slip),
        "-y" => CoulombFriction(:sticky),
        "+y" => CoulombFriction(:slip),
    ]
    for rigidbody in input.RigidBody
        @assert rigidbody.control == true
        # for calculation of effective mass in collision
        rigidbody.density = Inf
        rigidbody.model.m = Inf
    end
    @assert isempty(input.BoundaryCondition.Dirichlet)
end

function initialize(input::Input)
    GridState = MarbleRun.gridstate_type(input, Val(2), Float64)
    PointState = MarbleRun.pointstate_type(input, Val(2), Float64)

    # General
    coordsystem = input.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = input.General.domain
    dx = input.General.grid_space
    g = input.General.gravity

    # SoilLayer
    soillayers = input.SoilLayer
    H = sum(layer -> layer.thickness, soillayers) # ground surface
    @assert ymin + H ≤ ymax

    # Advanced
    α = input.Advanced.contact_threshold_scale
    nptsincell = input.Advanced.npoints_in_cell
    random_pts_gen = input.Advanced.random_points_generation

    grid = Grid(coordsystem, xmin:dx:xmax, ymin:dx:ymax)
    pointstate = generate_pointstate((x,y) -> y < ymin + H, PointState, grid; n=nptsincell, random=random_pts_gen)
    gridstate = generate_gridstate(GridState, grid)
    rigidbody = only(input.RigidBody).model

    bottom = ymin
    for i in length(soillayers):-1:1 # from low to high
        layer = soillayers[i]
        for p in 1:length(pointstate)
            y = pointstate.x[p][2]
            if bottom ≤ y ≤ bottom + layer.thickness
                pointstate.matindex[p] = i
            end
        end
        bottom += layer.thickness
    end

    # initialize variables of points
    for p in 1:length(pointstate)
        layerindex = pointstate.matindex[p]
        σ_y = 0.0
        for layer in soillayers[begin:layerindex-1]
            ρ₀ = layer.density
            σ_y += -ρ₀ * g * layer.thickness
        end
        h = sum(layer -> layer.thickness, soillayers[layerindex:end])
        layer = soillayers[layerindex]
        y = pointstate.x[p][2]
        ρ0 = layer.density
        K0 = layer.K0
        σ_y += -ρ0 * g * (h - y)
        σ_x = K0 * σ_y
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
        pointstate.m[p] = ρ0 * pointstate.V[p]
    end
    @. pointstate.x0 = pointstate.x
    @. pointstate.b = Vec(0.0, -g)

    if only(input.RigidBody).reset_position
        y0 = minimum(x -> x[2], coordinates(rigidbody))
        δ = sqrt(eps(Float64))
        translate!(rigidbody, Vec(0.0, (ymin - y0) + H + (α-1)*(dx/nptsincell)/2 + δ))
    else
        MarbleRun.remove_invalid_pointstate!(pointstate, input)
    end

    t = 0.0
    t, grid, gridstate, pointstate, rigidbody, deepcopy(rigidbody)
end

function main(input::Input, phase::Input_Phase, t, grid::Grid, gridstate::AbstractArray, pointstate::AbstractVector, rigidbody, rigidbody0)

    # General/Output
    dx = input.General.grid_space
    t_start = t
    t_stop = phase.time_stop
    t_step = input.Output.time_interval

    # SoilLayer
    matmodels = map(x -> x.model, input.SoilLayer)

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
    if input.Output.history
        outputs["history_file"] = joinpath(outdir, "history.csv")
        open(outputs["history_file"], "w") do io
            write(io, "depth,force\n")
        end
    end
    if input.Output.snapshots || input.Output.snapshot_first || input.Output.snapshot_last
        mkpath(joinpath(outdir, "snapshots"))
    end
    if isdefined(input.Injection, :main_output)
        Base.invokelatest(input.Injection.main_output_initialize, (;
            input,
            t,
            grid,
            gridstate,
            pointstate,
            rigidbody,
            rigidbody0,
        ))
    end

    ##################
    # Run simulation #
    ##################

    space = MPSpace(input.General.interpolation, grid, pointstate.x)
    logger = Logger(t_start, t_stop, t_step; input.General.showprogress)
    update!(logger, t)
    writeoutput(outputs, input, grid, gridstate, pointstate, rigidbody, rigidbody0, t, logindex(logger))

    if input.Output.snapshot_first
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot_first"),
            (; t, grid, gridstate, pointstate, rigidbody, rigidbody0)
        )
    end

    try
        while !isfinised(logger, t)
            dt = phase.CFL * MarbleRun.safe_minimum(LazyRows(pointstate)) do pt
                MarbleRun.timestep(matmodels[pt.matindex], pt, dx)
            end
            MarbleRun.advancestep!(grid, gridstate, pointstate, [rigidbody], space, dt, input, phase)

            if input.Output.quickview
                update!(logger, t += dt; print = MarbleRun.quickview_sparsity_pattern(@. !iszero($TransparentArray(gridstate.m))))
            else
                update!(logger, t += dt)
            end

            if islogpoint(logger)
                if input.Advanced.reorder_pointstate
                    Marble.reorder_pointstate!(pointstate, space)
                end
                writeoutput(outputs, input, grid, gridstate, pointstate, rigidbody, rigidbody0, t, logindex(logger))
            end
        end
    catch e
        writeoutput(outputs, input, grid, gridstate, pointstate, rigidbody, rigidbody0, t, "error")
        rethrow()
    end

    if input.Output.snapshot_last
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot_last"),
            (; t, grid, gridstate, pointstate, rigidbody, rigidbody0)
        )
    end

    t
end

function writeoutput(
        outputs::Dict{String, Any},
        input::Input,
        grid::Grid,
        gridstate::AbstractArray,
        pointstate::AbstractVector,
        rigidbody::GeometricObject,
        rigidbody0::GeometricObject,
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
                vtk_grid(vtm, rigidbody)
                if input.Paraview.GridState !== nothing
                    vtk_grid(vtm, grid; compress) do vtk
                        MarbleRun.writevtk(vtk, input.Paraview.GridState, gridstate)
                    end
                end
                pvd[t] = vtm
            end
        end
    end

    if input.Output.history
        history_file = outputs["history_file"]
        open(history_file, "a") do io
            (_, _), (ymin, _) = input.General.domain
            H = ymin + sum(layer -> layer.thickness, input.SoilLayer) # ground surface
            tip = minimum(x -> x[2], coordinates(rigidbody))
            depth = H - tip
            force = -sum(gridstate.fc)[2]
            if input.General.coordinate_system isa Axisymmetric
                force *= 2π
            end
            write(io, "$depth,$force\n")
        end
    end

    if input.Output.snapshots
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot$output_index"),
            (; t, grid, gridstate, pointstate, rigidbody, rigidbody0)
        )
    end

    if isdefined(input.Injection, :main_output)
        Base.invokelatest(input.Injection.main_output, (;
            input,
            grid,
            gridstate,
            pointstate,
            rigidbody,
            rigidbody0,
            t,
            output_index,
        ))
    end
end

end

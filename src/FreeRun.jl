module FreeRun

using PoingrSimulator
using PoingrSimulator: Input, Input_Phase
using Poingr
using GeometricObjects

using Serialization

struct NodeState
    m::Float64
    m′::Float64
    v::Vec{2, Float64}
    v_n::Vec{2, Float64}
    m_contacted::Float64
    vᵣ::Vec{2, Float64}
    fc::Vec{2, Float64}
    d::Vec{2, Float64}
    μ::Float64
end

struct PointState
    m::Float64
    V::Float64
    x::Vec{2, Float64}
    v::Vec{2, Float64}
    b::Vec{2, Float64}
    fc::Vec{2, Float64}
    σ::SymmetricSecondOrderTensor{3, Float64, 6}
    ϵ::SymmetricSecondOrderTensor{3, Float64, 6}
    ∇v::SecondOrderTensor{3, Float64, 9}
    C::Mat{2, 3, Float64, 6}
    r::Vec{2, Float64}
    index::Int
    matindex::Int
end

function preprocess_input!(input::Input)
    for mat in input.Material
        if hasproperty(mat, :friction_with_rigidbodies)
            @assert length(mat.friction_with_rigidbodies) == length(input.RigidBody)
        end
    end
end

function initialize(input::Input)
    coordinate_system = input.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = input.General.domain
    dx = input.General.grid_space
    g = input.General.gravity

    materials = input.Material
    rigidbodies = map(x -> x.model, input.RigidBody)

    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate(PointState, grid, input) do pointstate, matindex
        mat = materials[matindex]
        ρ0 = mat.density
        @. pointstate.m = ρ0 * pointstate.V
        @. pointstate.b = Vec(0.0, -g)
        @. pointstate.matindex = matindex
        PoingrSimulator.initialize_stress!(pointstate.σ, mat, g)
    end
    t = 0.0

    t, grid, pointstate, rigidbodies
end

function main(input::Input, phase::Input_Phase, t, grid, pointstate, rigidbodies)

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
    if input.Output.snapshots
        mkpath(joinpath(outdir, "snapshots"))
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
    logger = Logger(t_start, t_stop, t_step; input.General.show_progress)
    update!(logger, t)
    writeoutput(outputs, input, grid, pointstate, rigidbodies, t, logindex(logger))

    while !isfinised(logger, t)
        dt = phase.CFL * PoingrSimulator.safe_minimum(pointstate) do pt
            PoingrSimulator.timestep(matmodels[pt.matindex], pt, dx)
        end
        PoingrSimulator.advancestep!(grid, pointstate, rigidbodies, cache, dt, input, phase)
        update!(logger, t += dt)
        if islogpoint(logger)
            if input.Advanced.reorder_pointstate
                Poingr.reorder_pointstate!(pointstate, cache)
            end
            writeoutput(outputs, input, grid, pointstate, rigidbodies, t, logindex(logger))
        end
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
        output_index::Int,
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
                        vtk["nodal friction"] = vec(grid.state.μ)
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

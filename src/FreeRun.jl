module FreeRun

using PoingrSimulator
using PoingrSimulator: Input, getoftype
using Poingr
using GeometricObjects

using JLD2

struct NodeState
    m::Float64
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
    μ::Float64
    index::Int
    matindex::Int
end

function preprocess_input!(dict::Dict)
    get!(dict, "RigidBody", Ref(GeometricObject[]))
end

function initialize(INPUT::Input{:Root})
    coordinate_system = INPUT.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity

    materials = INPUT.Material
    rigidbodies = map(PoingrSimulator.create_rigidbody, INPUT.RigidBody)

    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate(PointState, grid, INPUT) do pointstate, matindex
        mat = materials[matindex]
        ρ0 = mat.density
        @. pointstate.m = ρ0 * pointstate.V
        @. pointstate.b = Vec(0.0, -g)
        @. pointstate.matindex = matindex
        PoingrSimulator.initialize_stress!(pointstate.σ, mat, g)
        if !isempty(rigidbodies)
            @. pointstate.μ = mat.friction_with_rigidbody
        end
    end
    t = 0.0

    grid, pointstate, rigidbodies, t
end

function main(INPUT::Input{:Root}, grid, pointstate, rigidbodies, t)

    #############
    # Constants #
    #############

    # General
    dx = INPUT.General.grid_space
    total_time = INPUT.General.total_time

    # Material models
    matmodels = map(PoingrSimulator.create_materialmodel, INPUT.Material)

    ################
    # Output files #
    ################

    output_dir = INPUT.Output.directory
    outputs = Dict{String, Any}()
    if INPUT.Output.snapshots
        outputs["snapshots_file"] = joinpath(output_dir, "snapshots.jld2")
        jldopen(identity, outputs["snapshots_file"], "w"; compress = true)
    end
    if INPUT.Output.paraview
        mkpath(joinpath(output_dir, "paraview"))
        outputs["paraview file"] = joinpath(output_dir, "paraview", "output")
        paraview_collection(vtk_save, outputs["paraview file"])
    end

    println("Particles: ", length(pointstate))

    ##################
    # Run simulation #
    ##################

    cache = MPCache(grid, pointstate.x)
    logger = Logger(0.0:INPUT.Output.interval:total_time; INPUT.General.show_progress)
    update!(logger, t)
    writeoutput(outputs, grid, pointstate, rigidbodies, logindex(logger), t, INPUT)

    while !isfinised(logger, t)
        dt = INPUT.Advanced.CFL * minimum(pointstate) do pt
            PoingrSimulator.timestep(matmodels[pt.matindex], pt, dx)
        end
        PoingrSimulator.advancestep!(grid, pointstate, rigidbodies, cache, INPUT, dt)
        update!(logger, t += dt)
        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)
            writeoutput(outputs, grid, pointstate, rigidbodies, logindex(logger), t, INPUT)
        end
    end
end

function writeoutput(
        outputs::Dict{String, Any},
        grid::Grid,
        pointstate::AbstractVector,
        rigidbodies::Vector{<: GeometricObject},
        output_index::Int,
        t::Real,
        INPUT::Input{:Root},
    )
    if INPUT.Output.paraview
        compress = true
        paraview_file = outputs["paraview file"]
        paraview_collection(paraview_file, append = true) do pvd
            vtk_multiblock(string(paraview_file, output_index)) do vtm
                vtk_points(vtm, pointstate.x; compress) do vtk
                    PoingrSimulator.write_vtk_points(vtk, pointstate)
                end
                for rigidbody in rigidbodies
                    vtk_grid(vtm, rigidbody)
                end
                if INPUT.Output.paraview_grid
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

    if INPUT.Output.snapshots
        jldopen(outputs["snapshots_file"], "a"; compress = true) do file
            file[string(output_index)] = (; grid, pointstate, rigidbodies, t)
        end
    end

    if isdefined(INPUT.Injection, :main_output)
        args = (;
            grid,
            pointstate,
            rigidbodies,
            INPUT,
            t,
            output_index,
        )
        INPUT.Injection.main_output(args)
    end
end

end

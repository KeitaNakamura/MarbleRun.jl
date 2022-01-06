module FreeRun

using PoingrSimulator
using PoingrSimulator: Input, getoftype
using Poingr
using GeometricObjects

using DelimitedFiles
using Serialization

struct NodeState
    m::Float64
    v::Vec{2, Float64}
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

function main(proj_dir::AbstractString, INPUT::Input{:Root}, Injection::Module)

    # General
    coordinate_system = INPUT.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity
    total_time = INPUT.General.total_time

    # Material
    materials = INPUT.Material

    # RigidBody
    rigidbodies = map(PoingrSimulator.create_rigidbody, get(INPUT, :RigidBody, GeometricObject[])::Vector)

    # Advanced
    α = getoftype(INPUT.Advanced, :contact_threshold_scale, 1.0)

    ##################
    # Initialization #
    ##################

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
    cache = MPCache(grid, pointstate.x)

    # constitutive models
    matmodels = map(materials) do mat
        PoingrSimulator.create_materialmodel(mat.type, mat, coordinate_system)
    end

    # boundary contacts
    boundary_contacts = PoingrSimulator.create_boundary_contacts(INPUT.BoundaryCondition)

    ################
    # Output files #
    ################

    outputs = Dict{String, Any}()
    # output directory
    output_dir = joinpath(proj_dir, INPUT.Output.folder_name)
    outputs["output directory"] = output_dir
    if INPUT.Output.serialize
        mkpath(joinpath(output_dir, "serialize"))
    end
    if INPUT.Output.paraview
        mkpath(joinpath(output_dir, "paraview"))
        outputs["paraview file"] = joinpath(output_dir, "paraview", "output")
        paraview_collection(vtk_save, outputs["paraview file"])
    end

    println("Particles: ", length(pointstate))

    t = 0.0
    logger = Logger(0.0:INPUT.Output.interval:total_time; INPUT.General.show_progress)
    update!(logger, t)
    writeoutput(outputs, grid, pointstate, rigidbodies, logindex(logger), t, INPUT, Injection)
    while !isfinised(logger, t)
        dt = INPUT.Advanced.CFL * minimum(pointstate) do pt
            PoingrSimulator.timestep(matmodels[pt.matindex], pt, dx)
        end

        update!(cache, grid, pointstate)
        PoingrSimulator.P2G!(grid, pointstate, cache, dt)
        masks = map(rigidbodies) do rigidbody
            PoingrSimulator.P2G_contact!(grid, pointstate, cache, dt, rigidbody, α, INPUT.Advanced.contact_penalty_parameter)
        end
        for bd in eachboundary(grid)
            @inbounds grid.state.v[bd.I] = boundary_velocity(grid.state.v[bd.I], bd.n, boundary_contacts)
        end
        PoingrSimulator.G2P!(pointstate, grid, cache, matmodels, materials, dt) # `materials` are for densities
        for (rigidbody, mask) in zip(rigidbodies, masks)
            PoingrSimulator.G2P_contact!(pointstate, grid, cache, rigidbody, mask)
        end

        for (i, (rigidbody, mask)) in enumerate(zip(rigidbodies, masks))
            input = INPUT.RigidBody[i]
            b = getoftype(input, :body_force, zero(Vec{2}))
            if getoftype(input, :control, false)
                GeometricObjects.update!(rigidbody, b, zero(Vec{3}), dt)
            else
                inds = findall(mask)
                Fc, Mc = GeometricObjects.compute_force_moment(rigidbody, view(pointstate.fc, inds), view(pointstate.x, inds))
                Fc += rigidbody.m * Vec(0,-g) + b
                GeometricObjects.update!(rigidbody, Fc, Mc, dt)
            end
        end

        update!(logger, t += dt)

        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)
            writeoutput(outputs, grid, pointstate, rigidbodies, logindex(logger), t, INPUT, Injection)
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
        Injection::Module,
    )
    if INPUT.Output.paraview
        paraview_file = outputs["paraview file"]
        paraview_collection(paraview_file, append = true) do pvd
            vtk_multiblock(string(paraview_file, output_index)) do vtm
                vtk_points(vtm, pointstate.x) do vtk
                    PoingrSimulator.write_vtk_points(vtk, pointstate)
                end
                for rigidbody in rigidbodies
                    vtk_grid(vtm, rigidbody)
                end
                if INPUT.Output.paraview_grid
                    vtk_grid(vtm, grid) do vtk
                        vtk["nodal contact force"] = vec(grid.state.fc)
                        vtk["nodal contact distance"] = vec(grid.state.d)
                        vtk["nodal friction"] = vec(grid.state.μ)
                    end
                end
                pvd[t] = vtm
            end
        end
    end

    if INPUT.Output.serialize
        serialize(joinpath(outputs["output directory"], "serialize", string("save", output_index)),
                  (; pointstate, grid, rigidbodies, t))
    end

    if isdefined(Injection, :main_output)
        args = (;
            grid,
            pointstate,
            rigidbodies,
            INPUT,
            t,
            output_index,
            output_dir = outputs["output directory"],
        )
        Injection.main_output(args)
    end
end

function boundary_velocity(v::Vec{2}, n::Vec{2}, boundary_contacts)
    n == Vec(-1,  0) && (side = :left)
    n == Vec( 1,  0) && (side = :right)
    n == Vec( 0, -1) && (side = :bottom)
    n == Vec( 0,  1) && (side = :top)
    contact = boundary_contacts[side]
    v + contact(v, n)
end

end

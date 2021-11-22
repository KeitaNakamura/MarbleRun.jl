module PenetrateIntoGround

using PoingrSimulator
using Poingr
using GeometricObjects

using DelimitedFiles
using Serialization
using TOML
using Dates

struct NodeState
    f::Vec{2, Float64}
    fc::Vec{2, Float64}
    fc_nor::Vec{2, Float64}
    w::Float64
    m::Float64
    v::Vec{2, Float64}
    vᵣ::Vec{2, Float64}
    w_rigidbody::Float64
    μ::Float64
end

struct PointState
    m::Float64
    V::Float64
    x::Vec{2, Float64}
    v::Vec{2, Float64}
    b::Vec{2, Float64}
    σ::SymmetricSecondOrderTensor{3, Float64, 6}
    ϵ::SymmetricSecondOrderTensor{3, Float64, 6}
    ∇v::SecondOrderTensor{3, Float64, 9}
    C::Mat{2, 3, Float64, 6}
    r::Vec{2, Float64}
    μ::Vector{Float64}
    index::Int
    matindex::Int
end

function julia_main()::Cint
    if isempty(ARGS)
        inputtoml = "input.toml"
    else
        inputtoml = ARGS[1]
    end
    try
        main(inputtoml)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function main(inputtoml::AbstractString)
    main(splitdir(inputtoml)[1],
         PoingrSimulator.parseinput(inputtoml))
end

function main(proj_dir::AbstractString, INPUT::NamedTuple)

    # General
    coordinate_system = INPUT.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity
    total_time = INPUT.General.total_time

    # SoilLayer
    soillayers = INPUT.SoilLayer # reorder layers from low to high
    H = sum(layer -> layer.thickness, soillayers)
    @assert H ≤ ymax

    # RigidBody
    rigidbody = Polygon(Vec{2}.(INPUT.RigidBody.coordinates))
    v_rigidbody = Vec(0.0, -INPUT.RigidBody.velocity)

    # Advanced
    α = INPUT.Advanced.contact_threshold_scale
    nptsincell = INPUT.Advanced.npoints_in_cell


    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate((x,y) -> y < H, PointState, grid; n = nptsincell)
    cache = MPCache(grid, pointstate.x)

    ##################
    # Initialization #
    ##################

    # constitutive models of layers
    matmodels = map(soillayers) do layer
        E = layer.youngs_modulus
        ν = layer.poissons_ratio
        c = layer.cohesion
        ϕ = layer.friction_angle
        ψ = layer.dilatancy_angle
        if layer.tension_cutoff === false
            tension_cutoff = Inf
        else
            tension_cutoff = layer.tension_cutoff
        end
        elastic = LinearElastic(; E, ν)
        if coordinate_system == "plane_strain"
            DruckerPrager(elastic, :plane_strain; c, ϕ, ψ, tension_cutoff)
        else
            DruckerPrager(elastic, :circumscribed; c, ϕ, ψ, tension_cutoff)
        end
    end

    # layer indices of points
    bottom = ymin
    for i in length(soillayers):-1:1 # from low to high
        layer = soillayers[i]
        Threads.@threads for p in 1:length(pointstate)
            y = pointstate.x[p][2]
            if bottom ≤ y ≤ bottom + layer.thickness
                pointstate.matindex[p] = i
            end
        end
        bottom += layer.thickness
    end

    # initialize variables of points
    Threads.@threads for p in 1:length(pointstate)
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
        ν = layer.poissons_ratio
        σ_y += -ρ0 * g * (h - y)
        σ_x = σ_y * ν / (1 - ν)
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
        pointstate.m[p] = ρ0 * pointstate.V[p]
        if layer.friction_with_rigidbody isa Vector
            pointstate.μ[p] = layer.friction_with_rigidbody
        else
            pointstate.μ[p] = [layer.friction_with_rigidbody]
        end
    end
    @. pointstate.b = Vec(0.0, -g)
    Poingr.reorder_pointstate!(pointstate, cache)

    translate!(rigidbody, Vec(0.0, H + (α-1)*(dx/nptsincell)/2))
    rigidbody_center_0 = centroid(rigidbody)

    ################
    # Output files #
    ################

    outputs = Dict{String, Any}()
    output_dir = joinpath(proj_dir, INPUT.Output.folder_name)
    mkpath(output_dir)
    outputs["output directory"] = output_dir
    ## paraview
    outputs["paraview file"] = joinpath(output_dir, "out")
    paraview_collection(vtk_save, outputs["paraview file"])
    ## history
    outputs["history file"] = joinpath(output_dir, "history.csv")
    open(outputs["history file"], "w") do io
        writedlm(io, ["disp" "force"], ',')
    end

    println("Start: ", now())
    println("Particles: ", length(pointstate))

    t = 0.0
    logger = Logger(0.0:INPUT.Output.interval:total_time; progress = INPUT.General.show_progress)
    while !isfinised(logger, t)
        dt = minimum(pointstate) do p
            ρ = p.m / p.V
            elastic = matmodels[p.matindex].elastic
            vc = matcalc(Val(:sound_speed), elastic.K, elastic.G, ρ)
            INPUT.Advanced.CFL * dx / vc
        end

        update!(cache, grid, pointstate)
        PoingrSimulator.P2G!(grid, pointstate, cache, dt)
        PoingrSimulator.P2G_contact!(grid, pointstate, cache, dt, rigidbody, v_rigidbody, α, INPUT.Advanced.contact_penalty_parameter)
        for bd in eachboundary(grid)
            @inbounds grid.state.v[bd.I] = boundary_velocity(grid.state.v[bd.I], bd.n)
        end
        PoingrSimulator.G2P!(pointstate, grid, cache, matmodels, dt)

        translate!(rigidbody, v_rigidbody * dt)
        update!(logger, t += dt)

        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)
            writeoutput(outputs, grid, pointstate, rigidbody, logger, rigidbody_center_0, t, INPUT)
        end
    end
end

function writeoutput(outputs::Dict{String, Any}, grid::Grid, pointstate::AbstractVector, rigidbody::Polygon, logger, rigidbody_center_0::Vec, t::Real, INPUT::NamedTuple)
    output_dir = outputs["output directory"]
    paraview_file = outputs["paraview file"]
    history_file = outputs["history file"]

    paraview_collection(paraview_file, append = true) do pvd
        vtk_multiblock(string(paraview_file, logindex(logger))) do vtm
            vtk_points(vtm, pointstate.x) do vtk
                PoingrSimulator.write_vtk_points(vtk, pointstate)
            end
            vtk_grid(vtm, rigidbody)
            if INPUT.Output.grid_in_paraview
                vtk_grid(vtm, grid) do vtk
                    vtk["nodal force"] = vec(grid.state.f)
                    vtk["nodal contact force"] = vec(grid.state.fc)
                    vtk["nodal contact force (normal)"] = vec(grid.state.fc_nor)
                    vtk["nodal friction"] = vec(grid.state.μ)
                end
            end
            pvd[t] = vtm
        end
    end

    open(history_file, "a") do io
        disp = abs(centroid(rigidbody)[2] - rigidbody_center_0[2])
        force = -sum(grid.state.fc)[2]
        if INPUT.General.coordinate_system == "axisymmetric"
            force *= 2π
        end
        writedlm(io, [disp force], ',')
    end

    serialize(joinpath(output_dir, string("save", logindex(logger))),
              (; pointstate, grid, rigidbody))
end

function boundary_velocity(v::Vec, n::Vec)
    if n == Vec(0, -1) # bottom
        v + Contact(:sticky)(v, n)
    else
        v + Contact(:slip)(v, n)
    end
end

end

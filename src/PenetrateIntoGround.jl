module PenetrateIntoGround

using PoingrSimulator
using PoingrSimulator: Input
using Poingr
using GeometricObjects

using Serialization

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
    σ::SymmetricSecondOrderTensor{3, Float64, 6}
    ϵ::SymmetricSecondOrderTensor{3, Float64, 6}
    ∇v::SecondOrderTensor{3, Float64, 9}
    C::Mat{2, 3, Float64, 6}
    r::Vec{2, Float64}
    μ::Vector{Float64}
    index::Int
    matindex::Int
end

function preprocess_input!(dict::Dict)
    dict["Material"] = dict["SoilLayer"]
    for mat in dict["Material"]
        mat["type"] = DruckerPrager
    end
    dict["BoundaryCondition"] = Dict{String, Any}("bottom" => Inf)
    for rigidbody in dict["RigidBody"]
        rigidbody["density"] = Inf # for calculation of effective mass in collision
    end
end

function initialize(INPUT::Input{:Root})
    # General
    coordinate_system = INPUT.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity

    # SoilLayer
    soillayers = INPUT.SoilLayer
    H = sum(layer -> layer.thickness, soillayers) # ground surface
    @assert H ≤ ymax

    # Advanced
    α = INPUT.Advanced.contact_threshold_scale
    nptsincell = INPUT.Advanced.npoints_in_cell

    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate((x,y) -> y < H, PointState, grid; n = nptsincell)
    rigidbody = PoingrSimulator.create_rigidbody(only(INPUT.RigidBody))

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
        if layer.friction_with_rigidbody isa Tuple
            pointstate.μ[p] = collect(layer.friction_with_rigidbody)
        else
            pointstate.μ[p] = [layer.friction_with_rigidbody]
        end
    end
    @. pointstate.b = Vec(0.0, -g)

    translate!(rigidbody, Vec(0.0, H + (α-1)*(dx/nptsincell)/2))
    t = 0.0

    grid, pointstate, rigidbody, deepcopy(rigidbody), t
end

function main(INPUT::Input{:Root}, grid, pointstate, rigidbody, rigidbody0, t)

    println("Particles: ", length(pointstate))

    # General/Output
    dx = INPUT.General.grid_space
    t_stop = INPUT.General.total_time
    t_step = INPUT.Output.interval
    t_start = (0.0:t_step:t_stop)[searchsortedlast(0.0:t_step:t_stop, t)]

    # SoilLayer
    matmodels = map(PoingrSimulator.create_materialmodel, INPUT.SoilLayer)

    ################
    # Output files #
    ################

    output_dir = INPUT.Output.directory
    outputs = Dict{String, Any}()
    if INPUT.Output.paraview
        mkpath(joinpath(output_dir, "paraview"))
        outputs["paraview_file"] = joinpath(output_dir, "paraview", "output")
        paraview_collection(vtk_save, outputs["paraview_file"])
    end
    if INPUT.Output.history
        outputs["history_file"] = joinpath(output_dir, "history.csv")
        open(outputs["history_file"], "w") do io
            write(io, "disp,force\n")
        end
    end
    if INPUT.Output.snapshots
        mkpath(joinpath(output_dir, "snapshots"))
    end
    if isdefined(INPUT.Injection, :main_output)
        INPUT.Injection.main_output_initialize((;
            INPUT,
            grid,
            pointstate,
            rigidbody,
            rigidbody0,
            t,
        ))
    end

    ##################
    # Run simulation #
    ##################

    cache = MPCache(grid, pointstate.x)
    logger = Logger(t_start:t_step:t_stop; INPUT.General.show_progress)
    update!(logger, t)
    writeoutput(outputs, INPUT, grid, pointstate, rigidbody, rigidbody0, t, logindex(logger))

    while !isfinised(logger, t)
        dt = INPUT.Advanced.CFL * minimum(pointstate) do pt
            PoingrSimulator.timestep(matmodels[pt.matindex], pt, dx)
        end
        PoingrSimulator.advancestep!(grid, pointstate, [rigidbody], cache, INPUT, dt)
        update!(logger, t += dt)
        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)
            writeoutput(outputs, INPUT, grid, pointstate, rigidbody, rigidbody0, t, logindex(logger))
        end
    end
end

function writeoutput(
        outputs::Dict{String, Any},
        INPUT::Input{:Root},
        grid::Grid,
        pointstate::AbstractVector,
        rigidbody::GeometricObject,
        rigidbody0::GeometricObject,
        t::Real,
        output_index::Int,
    )
    if INPUT.Output.paraview
        compress = true
        paraview_file = outputs["paraview_file"]
        paraview_collection(paraview_file, append = true) do pvd
            vtk_multiblock(string(paraview_file, output_index)) do vtm
                vtk_points(vtm, pointstate.x; compress) do vtk
                    PoingrSimulator.write_vtk_points(vtk, pointstate)
                end
                vtk_grid(vtm, rigidbody)
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

    if INPUT.Output.history
        history_file = outputs["history_file"]
        open(history_file, "a") do io
            disp = abs(centroid(rigidbody)[2] - centroid(rigidbody0)[2])
            force = -sum(grid.state.fc)[2]
            if INPUT.General.coordinate_system isa Axisymmetric
                force *= 2π
            end
            write(io, "$disp,$force\n")
        end
    end

    if INPUT.Output.snapshots
        serialize(
            joinpath(INPUT.Output.directory, "snapshots", "snapshot$output_index"),
            (; grid, pointstate, rigidbody, rigidbody0, t)
        )
    end

    if isdefined(INPUT.Injection, :main_output)
        INPUT.Injection.main_output((;
            INPUT,
            grid,
            pointstate,
            rigidbody,
            rigidbody0,
            t,
            output_index,
        ))
    end
end

end

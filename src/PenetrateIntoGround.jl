module PenetrateIntoGround

using PoingrSimulator
using PoingrSimulator: Input
using Poingr
using GeometricObjects

using DelimitedFiles
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

function main(proj_dir::AbstractString, INPUT::Input{:Root}, Injection::Module)

    # General
    coordinate_system = INPUT.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity
    total_time = INPUT.General.total_time

    # SoilLayer
    soillayers = INPUT.SoilLayer
    H = sum(layer -> layer.thickness, soillayers)
    @assert H ≤ ymax

    # RigidBody
    rigidbody = PoingrSimulator.create_rigidbody(only(INPUT.RigidBody))

    # Advanced
    α = INPUT.Advanced.contact_threshold_scale
    nptsincell = INPUT.Advanced.npoints_in_cell


    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate((x,y) -> y < H, PointState, grid; n = nptsincell)
    cache = MPCache(grid, pointstate.x)

    ##################
    # Initialization #
    ##################

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
        if layer.friction_with_rigidbody isa Tuple
            pointstate.μ[p] = collect(layer.friction_with_rigidbody)
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
    # output directory
    output_dir = joinpath(proj_dir, INPUT.Output.folder_name)
    outputs["output directory"] = output_dir
    if INPUT.Output.serialize
        outputs["serialized_data_file"] = joinpath(outputs["output directory"], "serialized_data.jld2")
        jldopen(identity, outputs["serialized_data_file"], "w"; compress = true)
    end
    if INPUT.Output.paraview
        mkpath(joinpath(output_dir, "paraview"))
        outputs["paraview file"] = joinpath(output_dir, "paraview", "output")
        paraview_collection(vtk_save, outputs["paraview file"])
    end
    if INPUT.Output.history
        outputs["history file"] = joinpath(output_dir, "history.csv")
        open(outputs["history file"], "w") do io
            writedlm(io, ["disp" "force"], ',')
        end
    end

    println("Particles: ", length(pointstate))

    t = 0.0
    logger = Logger(0.0:INPUT.Output.interval:total_time; INPUT.General.show_progress)
    update!(logger, t)
    writeoutput(outputs, grid, pointstate, rigidbody, logindex(logger), rigidbody_center_0, t, INPUT, Injection)
    while !isfinised(logger, t)
        dt = PoingrSimulator.advancestep!(grid, pointstate, [rigidbody], cache, INPUT)
        update!(logger, t += dt)
        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)
            writeoutput(outputs, grid, pointstate, rigidbody, logindex(logger), rigidbody_center_0, t, INPUT, Injection)
        end
    end
end

function writeoutput(
        outputs::Dict{String, Any},
        grid::Grid,
        pointstate::AbstractVector,
        rigidbody::GeometricObject,
        output_index::Int,
        rigidbody_center_0::Vec,
        t::Real,
        INPUT::Input{:Root},
        Injection::Module,
    )
    if INPUT.Output.paraview
        compress = true
        paraview_file = outputs["paraview file"]
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
        history_file = outputs["history file"]
        open(history_file, "a") do io
            disp = abs(centroid(rigidbody)[2] - rigidbody_center_0[2])
            force = -sum(grid.state.fc)[2]
            if INPUT.General.coordinate_system isa Axisymmetric
                force *= 2π
            end
            writedlm(io, [disp force], ',')
        end
    end

    if INPUT.Output.serialize
        jldopen(outputs["serialized_data_file"], "a"; compress = true) do file
            file[string(output_index)] = (; pointstate, grid, rigidbody, t)
        end
    end

    if isdefined(Injection, :main_output)
        args = (;
            grid,
            pointstate,
            rigidbody,
            INPUT,
            t,
            output_index,
            output_dir = outputs["output directory"],
        )
        Injection.main_output(args)
    end
end

function boundary_velocity(v::Vec, n::Vec)
    if n == Vec(0, -1) # bottom
        v + Contact(:sticky)(v, n)
    else
        v + Contact(:slip)(v, n)
    end
end

end

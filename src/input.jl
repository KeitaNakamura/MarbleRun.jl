struct Input{name, Tup <: NamedTuple}
    fields::Tup
end
Input{name}(fields::NamedTuple) where {name} = Input{name, typeof(fields)}(fields)
Input{name}(; kwargs...) where {name} = Input{name}((; kwargs...))

Base.propertynames(input::Input) = propertynames(getfield(input, :fields))
Base.getproperty(input::Input, name::Symbol) = getproperty(getfield(input, :fields), name)

Base.length(input::Input) = length(getfield(input, :fields))
Base.getindex(input::Input, name) = getindex(getfield(input, :fields), name)
Base.iterate(input::Input) = iterate(getfield(input, :fields))
Base.iterate(input::Input, state) = iterate(getfield(input, :fields), state)

Base.get(input::Input, name::Symbol, default) = get(getfield(input, :fields), name, default)
Base.haskey(input::Input, name::Symbol) = haskey(getfield(input, :fields), name)
Base.keys(input::Input) = keys(getfield(input, :fields))
Base.values(input::Input) = values(getfield(input, :fields))
Base.merge(tup::NamedTuple, input::Input) = merge(tup, getfield(input, :fields))
Base.show(io::IO, input::Input{name}) where {name} = print(io, "Input{:$name}", getfield(input, :fields))

getoftype(input::Input, name::Symbol, default)::typeof(default) = oftype(default, get(getfield(input, :fields), name, default))
getoftype(input::Dict, name, default)::typeof(default) = oftype(default, get(input, name, default))

##############
# Input TOML #
##############

_parse_input(name, x) = x
_parse_input(name, x::Base.RefValue) = x[]
_parse_input(name, x::Vector) = first(x) isa Dict ? _parse_input.(name, x) : (x...,) # try vector => tuple except for table
_parse_input(name, x::Dict) = Input{name}(; (Symbol(key) => _parse_input(Symbol(key), value) for (key, value) in x)...)

function parse_input(x::Dict)
    for section in keys(x)
        preprocess! = Symbol(:preprocess_, section, :!)
        isdefined(@__MODULE__, preprocess!) && eval(preprocess!)(x[section])
    end
    if haskey(x["General"], "type")
        x["General"]["type"].preprocess_input!(x)
    end
    _parse_input(:Root, x)
end

parse_inputfile(path::AbstractString) = parse_input(TOML.parsefile(path))
parse_inputstring(str::AbstractString) = parse_input(TOML.parse(str))


# helper functions for preprocess
function eval_convert(::Type{T}, str::String)::T where {T}
    convert(T, eval(Meta.parse(str)))
end
function eval_convert(::Type{T}, val)::T where {T}
    convert(T, val)
end
function ifhaskey_eval_convert!(::Type{T}, dict::Dict, name::String) where {T}
    if haskey(dict, name)
        dict[name] = eval_convert(T, dict[name])
    end
end

###########
# General #
###########

function preprocess_General!(General::Dict)
    ifhaskey_eval_convert!(Module, General, "type")

    if haskey(General, "coordinate_system")
        coordinate_system = General["coordinate_system"]
        if coordinate_system == "plane_strain"
            General["coordinate_system"] = PlaneStrain()
        elseif coordinate_system == "axisymmetric"
            General["coordinate_system"] = Axisymmetric()
        else
            throw(ArgumentError("wrong `coordinate_system`, got \"$coordinate_system\", use \"plane_strain\" or \"axisymmetric\""))
        end
    end
end

#####################
# BoundaryCondition #
#####################

function preprocess_BoundaryCondition!(BoundaryCondition::Dict)
    for side in ("left", "right", "bottom", "top")
        coef = getoftype(BoundaryCondition, side, 0.0)
        BoundaryCondition[side] = Contact(:friction, coef) # friction can also handle sticky and slip
    end
end

#############
# SoilLayer #
#############

function preprocess_SoilLayer!(SoilLayer::Vector)
    for layer in SoilLayer
        if haskey(layer, "friction_with_rigidbodies")
            wrap(x)::Vector{Float64} = x isa Number ? [x] : x
            layer["friction_with_rigidbodies"] = map(wrap, layer["friction_with_rigidbodies"]) # handle friction coefficients with polygon object
        end
    end
end

############
# Material #
############

const InputMaterial = Union{Input{:Material}, Input{:SoilLayer}}

function preprocess_Material!(Material::Vector)
    for mat in Material
        ifhaskey_eval_convert!(Function,               mat, "region")
        ifhaskey_eval_convert!(Type{<: MaterialModel}, mat, "type")
        if haskey(mat, "friction_with_rigidbodies")
            wrap(x)::Vector{Float64} = x isa Number ? [x] : x
            mat["friction_with_rigidbodies"] = map(wrap, mat["friction_with_rigidbodies"]) # handle friction coefficients with polygon object
        end
    end
end

function create_materialmodel(mat::InputMaterial)
    create_materialmodel(mat.type, mat)
end

function create_materialmodel(::Type{DruckerPrager}, params::InputMaterial)
    mohr_coulomb_type = params.mohr_coulomb_type
    E = params.youngs_modulus
    ν = params.poissons_ratio
    c = params.cohesion
    ϕ = params.friction_angle
    ψ = params.dilatancy_angle
    tension_cutoff = params.tension_cutoff
    elastic = LinearElastic(; E, ν)
    DruckerPrager(elastic, mohr_coulomb_type; c, ϕ, ψ, tension_cutoff)
end

function create_materialmodel(::Type{NewtonianFluid}, params::InputMaterial)
    ρ0 = params.density
    P0 = params.pressure
    c = params.sound_of_speed
    μ = params.viscosity
    NewtonianFluid(; ρ0, P0, c, μ)
end

# This function is basically based on Material.Initialization
function initialize_stress!(σₚ::AbstractVector, material::Input{:Material}, g)
    Initialization = material.Initialization
    ρ0 = material.density
    if Initialization.type == "K0"
        ν = material.poissons_ratio
        K0 = Initialization.K0 == "auto" ? ν / (1 - ν) : Initialization.K0
        for p in eachindex(σₚ)
            σ_y = -ρ0 * g * Initialization.reference_height
            σ_x = K0 * σ_y
            σₚ[p] = (@Mat [σ_x 0.0 0.0
                           0.0 σ_y 0.0
                           0.0 0.0 σ_x]) |> symmetric
        end
    elseif Initialization.type == "uniform"
        for p in eachindex(σₚ)
            σₚ[p] = Initialization.mean_stress * one(σₚ[p])
        end
    else
        throw(ArgumentError("invalid initialization type, got $(Initialization.type)"))
    end
end

# Since initializing pointstate is dependent on types of simulation,
# `initialize!` phase should be given as argument.
# This `initialize!` function is called on each material to perform
# material-wise initialization.
# After initialization, these `pointstate`s will be concatenated.
# Then, if you have rigid bodies, the invalid points are removed.
function Poingr.generate_pointstate(initialize!::Function, ::Type{PointState}, grid::Grid, INPUT::Input{:Root}) where {PointState}
    # generate all pointstate first
    Material = INPUT.Material
    pointstates = map(1:length(Material)) do matindex
        material = Material[matindex]
        pointstate′ = generate_pointstate( # call method in `Poingr`
            material.region,
            PointState,
            grid;
            n = getoftype(INPUT.Advanced, :npoints_in_cell, 2),
        )
        initialize!(pointstate′, matindex)
        pointstate′
    end
    pointstate = first(pointstates)
    append!(pointstate, pointstates[2:end]...)

    # remove invalid pointstate
    α = getoftype(INPUT.Advanced, :contact_threshold_scale, 1.0)
    !isempty(INPUT.RigidBody) && deleteat!(
        pointstate,
        findall(eachindex(pointstate)) do p
            xₚ = pointstate.x[p]
            rₚ = pointstate.r[p]
            any(map(create_rigidbody, INPUT.RigidBody)) do rigidbody
                # remove pointstate which is in rigidbody or is in contact with rigidbody
                in(xₚ, rigidbody) || distance(rigidbody, xₚ, α * mean(rₚ)) !== nothing
            end
        end,
    )

    pointstate
end

#############
# RigidBody #
#############

function preprocess_RigidBody!(RigidBody::Vector)
    foreach(preprocess_RigidBody!, RigidBody)
end

function preprocess_RigidBody!(RigidBody::Dict)
    ifhaskey_eval_convert!(Type{<: Shape}, RigidBody, "type")
    if haskey(RigidBody, "control")
        if RigidBody["control"] === true
            # If the rigid body is controled, `density` should be `Inf`
            # to set `mass` into `Inf`.
            # This is necessary for computing effective mass.
            RigidBody["density"] = Inf # TODO: warning?
        end
    end
end

function create_rigidbody(RigidBody::Input{:RigidBody})
    create_rigidbody(RigidBody.type, RigidBody)
end

function create_rigidbody(::Type{Polygon}, params::Input{:RigidBody})
    rigidbody = GeometricObject(Polygon(Vec{2}.(params.coordinates)...))
    initialize_rigidbody!(rigidbody, params)
    rigidbody
end

function create_rigidbody(::Type{Circle}, params::Input{:RigidBody})
    rigidbody = GeometricObject(Circle(Vec(params.center), params.radius))
    initialize_rigidbody!(rigidbody, params)
    rigidbody
end

function initialize_rigidbody!(rigidbody::GeometricObject{dim, T}, params::Input{:RigidBody}) where {dim, T}
    rigidbody.m = area(rigidbody) * params.density # TODO: should use `volume`?
    rigidbody.v = getoftype(params, :velocity, zero(Vec{dim, T}))
    rigidbody.ω = getoftype(params, :angular_velocity, zero(Vec{3, T}))
    rigidbody
end

##########
# Output #
##########

function preprocess_Output!(Output::Dict)
end

############
# Advanced #
############

function preprocess_Advanced!(Advanced::Dict)
end

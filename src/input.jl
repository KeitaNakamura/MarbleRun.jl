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
Base.show(io::IO, input::Input{name}) where {name} = print(io, "Input{:$name}", getfield(input, :fields))


##############
# Input TOML #
##############

_parse_input(name, x) = x
_parse_input(name, x::Vector) = first(x) isa Dict ? _parse_input.(name, x) : (x...,) # try vector => tuple except for table
_parse_input(name, x::Dict) = Input{name}(; (Symbol(key) => _parse_input(Symbol(key), value) for (key, value) in x)...)

function parse_input(x::Dict)
    for section in keys(x)
        preprocess! = Symbol(:preprocess_, section, :!)
        eval(preprocess!)(x[section])
    end
    _parse_input(:Root, x)
end

parse_inputfile(path::AbstractString) = parse_input(TOML.parsefile(path))
parse_inputstring(str::AbstractString) = parse_input(TOML.parse(str))

###########
# General #
###########

function preprocess_General!(General::Dict)
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
end

function create_boundary_contacts(BoundaryCondition::Input{:BoundaryCondition})
    dict = Dict{Symbol, Contact}()
    for side in (:left, :right, :bottom, :top)
        if haskey(BoundaryCondition, side)
            coef = BoundaryCondition[side]
            coef = convert(Float64, coef isa AbstractString ? eval(Meta.parse(coef)) : coef)
            contact = Contact(:friction, coef)
        else
            contact = Contact(:slip)
        end
        dict[side] = contact
    end
    dict
end

#############
# SoilLayer #
#############

function preprocess_SoilLayer!(SoilLayer::Vector)
end

############
# Material #
############

const InputMaterial = Union{Input{:Material}, Input{:SoilLayer}}

function preprocess_Material!(Material::Vector)
    for mat in Material
        if haskey(mat, "region")
            mat["region"] = eval(Meta.parse(mat["region"])) # should be anonymous function
        end
        if haskey(mat, "type")
            mat["type"] = eval(Meta.parse(mat["type"]))
        end
    end
end

function create_materialmodel(mat::InputMaterial, coordinate_system)
    create_materialmodel(mat.type, mat, coordinate_system)
end

function create_materialmodel(::Type{DruckerPrager}, params::InputMaterial, coordinate_system)
    E = params.youngs_modulus
    ν = params.poissons_ratio
    c = params.cohesion
    ϕ = params.friction_angle
    ψ = params.dilatancy_angle
    tension_cutoff = params.tension_cutoff
    elastic = LinearElastic(; E, ν)
    if coordinate_system isa PlaneStrain
        DruckerPrager(elastic, :plane_strain; c, ϕ, ψ, tension_cutoff)
    else
        DruckerPrager(elastic, :circumscribed; c, ϕ, ψ, tension_cutoff)
    end
end

function create_materialmodel(::Type{NewtonianFluid}, params::InputMaterial, coordinate_system)
    ρ0 = params.density
    P0 = params.pressure
    c = params.sound_of_speed
    μ = params.viscosity
    NewtonianFluid(; ρ0, P0, c, μ)
end

# Material.Initialization

function initialize_stress!(σₚ::AbstractVector, material::Input{:Material}, g)
    Initialization = material.Initialization
    ρ0 = material.density
    if Initialization.type == "K0"
        for p in eachindex(σₚ)
            σ_y = -ρ0 * g * Initialization.reference_height
            σ_x = Initialization.K0 * σ_y
            σₚ[p] = (@Mat [σ_x 0.0 0.0
                           0.0 σ_y 0.0
                           0.0 0.0 σ_x]) |> symmetric
        end
    elseif Initialization.type == "uniform"
        for p in eachindex(σₚ)
            σₚ[p] = Initialization.mean_stress * one(σₚ[p])
        end
    else
        throw(ArgumentError("invalid initialization type, got $(condition.type)"))
    end
end

#############
# RigidBody #
#############

function preprocess_RigidBody!(RigidBody::Dict)
end

function preprocess_RigidBody!(RigidBody::Vector)
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

using Poingr
using Poingr: Interpolation
using GeometricObjects
using TOML

function parse_input(str::AbstractString; project = ".")
    input = convert_input(TOMLInput(TOML.parse(str)))
    input.project = project
    input.Output.directory = joinpath(input.project, input.Output.directory)
    input.General.type.preprocess_input!(input)
    input
end
parse_inputfile(path::AbstractString) = parse_input(read(path, String); project = dirname(path))

abstract type TOMLTable end

struct EvalString{T}
    content::T
end
Base.convert(::Type{EvalString{T}}, x::String) where {T} = EvalString{T}(eval(Meta.parse(x)))
convert_input(x::EvalString) = x.content

struct ToVec
    content::Vector{Float64}
end
Base.convert(::Type{ToVec}, x::Vector) = ToVec(x)
convert_input(x::ToVec) = Vec{length(x.content), Float64}(x.content)

struct ToTuple{T}
    content::Vector{T}
end
Base.convert(::Type{ToTuple{T}}, x::Vector) where {T} = ToTuple{T}(x)
convert_input(x::ToTuple{T}) where {T} = NTuple{length(x.content), T}(x.content)

struct SkipEntry{T}
    content::T
end
Base.convert(::Type{SkipEntry{T}}, x) where {T} = SkipEntry{T}(x)

_undefkey(s) = throw(UndefKeywordError(s))

###########
# General #
###########

Base.@kwdef mutable struct TOMLInput_General <: TOMLTable
    type              :: EvalString{Module}
    coordinate_system :: String
    domain            :: Vector{Vector{Float64}}
    grid_space        :: Float64
    gravity           :: Float64
    interpolation     :: EvalString{Interpolation} = "LinearWLS(QuadraticBSpline())"
    transfer          :: ToTuple{String}           = startswith(interpolation, "LinearWLS") ||
                                                     startswith(interpolation, "KernelCorrection") ?
                                                     ["affine", "PIC"] : ["normal", "FLIP"]
    show_progress     :: Bool                      = true
end

mutable struct Input_General{CoordSystem <: CoordinateSystem, Interp <: Interpolation}
    type              :: Module
    coordinate_system :: CoordSystem
    domain            :: Vector{Vector{Float64}}
    grid_space        :: Float64
    gravity           :: Float64
    interpolation     :: Interp
    transfer          :: Tuple{String, String}
    show_progress     :: Bool
end

function convert_input(input::TOMLInput_General, ::Val{:coordinate_system})
    input.coordinate_system == "plane_strain" && return PlaneStrain()
    input.coordinate_system == "axisymmetric" && return Axisymmetric()
    throw(ArgumentError("wrong `coordinate_system`, got \"$coordinate_system\", use \"plane_strain\" or \"axisymmetric\""))
end


#########
# Phase #
#########

Base.@kwdef mutable struct TOMLInput_Phase <: TOMLTable
    time_stop     :: Float64
    CFL           :: Float64 = 0.5
    restart       :: String  = ""
    update_motion :: Bool    = true
end

mutable struct Input_Phase
    time_stop     :: Float64
    CFL           :: Float64
    restart       :: String
    update_motion :: Bool
end

#####################
# BoundaryCondition #
#####################

# dirichlet
Base.@kwdef struct TOMLInput_BoundaryCondition_Dirichlet <: TOMLTable
    inbounds       :: EvalString{Function}
    velocity       :: ToVec
    output         :: Bool                 = true
end
mutable struct Input_BoundaryCondition_Dirichlet{dim}
    inbounds       :: Function
    velocity       :: Vec{dim, Float64}
    output         :: Bool
    displacement   :: Float64
    reaction_force :: Float64
    active_nodes   :: Array{Bool, dim}
end
function convert_input(input::TOMLInput_BoundaryCondition_Dirichlet)
    inbounds = convert_input(input.inbounds)
    velocity = convert_input(input.velocity)
    output = input.output
    dim = length(velocity)
    Input_BoundaryCondition_Dirichlet(inbounds, velocity, output, 0.0, 0.0, Array{Bool, dim}(undef, fill(0, dim)...))
end

Base.@kwdef mutable struct TOMLInput_BoundaryCondition <: TOMLTable
    top    :: Float64 = 0.0
    bottom :: Float64 = 0.0
    left   :: Float64 = 0.0
    right  :: Float64 = 0.0
    Dirichlet :: Vector{TOMLInput_BoundaryCondition_Dirichlet} = []
end

mutable struct Input_BoundaryCondition{Dirichlet}
    top       :: Contact
    bottom    :: Contact
    left      :: Contact
    right     :: Contact
    Dirichlet :: Vector{Dirichlet}
end

function convert_input(input::TOMLInput_BoundaryCondition, ::Val{name}) where {name}
    name != :Dirichlet && return Contact(:friction, getproperty(input, name))
    map(convert_input, input.Dirichlet)
end

##########
# Output #
##########

Base.@kwdef mutable struct TOMLInput_Output <: TOMLTable
    time_interval  :: Float64
    directory      :: String  = "output.tmp"
    snapshots      :: Bool    = false
    paraview       :: Bool    = true
    paraview_grid  :: Bool    = false
    copy_inputfile :: Bool    = true
    history        :: Bool    = true # only for `PenetrateIntoGround`
end

mutable struct Input_Output
    time_interval  :: Float64
    directory      :: String
    snapshots      :: Bool
    paraview       :: Bool
    paraview_grid  :: Bool
    copy_inputfile :: Bool
    history        :: Bool
end

############
# Material #
############

const VV{T} = Vector{Vector{T}}
wrap(x)::Vector{Float64} = x isa Number ? [x] : x

Base.@kwdef mutable struct TOMLInput_Material <: TOMLTable
    region                              :: EvalString{Function}
    density                             :: Float64
    model                               :: MaterialModel
    init
    friction_with_rigidbodies           :: VV{Float64}            = []
    friction_threshold_with_rigidbodies :: SkipEntry{VV{Float64}} = []
    function TOMLInput_Material(region, density, model, init, friction_with_rigidbodies, friction_threshold_with_rigidbodies)
        new(region, density, model, init, map(wrap, friction_with_rigidbodies), map(wrap, friction_threshold_with_rigidbodies))
    end
end

mutable struct Input_Material{Model <: MaterialModel, Init}
    region                    :: Function
    density                   :: Float64
    model                     :: Model
    init                      :: Init
    friction_with_rigidbodies :: VV{Vec{2, Float64}} # Vec{2}(μ, c)
end

function convert_input(input::TOMLInput_Material, ::Val{:friction_with_rigidbodies})::VV{Vec{2, Float64}}
    μ = input.friction_with_rigidbodies
    c = input.friction_threshold_with_rigidbodies.content
    c = isempty(c) ? [fill(0.0, length(x)) for x in μ] : c
    [[Vec(y) for y in zip(x...)] for x in zip(μ, c)]
end

#############
# SoilLayer #
#############

Base.@kwdef mutable struct TOMLInput_SoilLayer <: TOMLTable
    thickness                           :: Float64
    density                             :: Float64
    poissons_ratio                      :: Float64                = NaN
    K0                                  :: Float64                = isnan(poissons_ratio) ? _undefkey(:K0) : poissons_ratio / (1 - poissons_ratio)
    model                               :: MaterialModel
    friction_with_rigidbodies           :: VV{Float64}            = []
    friction_threshold_with_rigidbodies :: SkipEntry{VV{Float64}} = []
    function TOMLInput_SoilLayer(region, density, poissons_ratio, K0, model, friction_with_rigidbodies, friction_threshold_with_rigidbodies)
        new(region, density, poissons_ratio, K0, model, map(wrap, friction_with_rigidbodies), map(wrap, friction_threshold_with_rigidbodies))
    end
end

mutable struct Input_SoilLayer{Model <: MaterialModel}
    thickness                 :: Float64
    density                   :: Float64
    poissons_ratio            :: Float64
    K0                        :: Float64
    model                     :: Model
    friction_with_rigidbodies :: VV{Vec{2, Float64}} # Vec{2}(μ, c)
end

function convert_input(input::TOMLInput_SoilLayer, ::Val{:friction_with_rigidbodies})::VV{Vec{2, Float64}}
    μ = input.friction_with_rigidbodies
    c = input.friction_threshold_with_rigidbodies.content
    c = isempty(c) ? [fill(0.0, length(x)) for x in μ] : c
    [[Vec(y) for y in zip(x...)] for x in zip(μ, c)]
end

##############################
# (Material/SoilLayer).model #
##############################

# newtonian fluid

Base.@kwdef mutable struct TOMLInput_Material_model_NewtonianFluid <: TOMLTable
    density_ref    :: Float64
    pressure_ref   :: Float64
    speed_of_sound :: Float64
    viscosity      :: Float64
end

function Base.convert(::Type{MaterialModel}, model::TOMLInput_Material_model_NewtonianFluid)
    NewtonianFluid(;
        ρ0 = model.density_ref,
        P0 = model.pressure_ref,
        c = model.speed_of_sound,
        μ = model.viscosity,
    )
end

# Drucker-Prager model

Base.@kwdef mutable struct TOMLInput_Material_model_DruckerPrager <: TOMLTable
    mohr_coulomb_type :: String
    youngs_modulus    :: Float64
    poissons_ratio    :: Float64
    cohesion          :: Float64
    friction_angle    :: Float64
    dilatancy_angle   :: Float64
    tension_cutoff    :: Float64 = 0.0
end

const TOMLInput_SoilLayer_model_DruckerPrager = TOMLInput_Material_model_DruckerPrager

function Base.convert(::Type{MaterialModel}, model::TOMLInput_Material_model_DruckerPrager)
    DruckerPrager(
        LinearElastic(; E = model.youngs_modulus, ν = model.poissons_ratio),
        model.mohr_coulomb_type;
        c = model.cohesion,
        ϕ = model.friction_angle,
        ψ = model.dilatancy_angle,
        tension_cutoff = model.tension_cutoff,
    )
end

#############################
# (Material/SoilLayer).init #
#############################

Base.@kwdef mutable struct TOMLInput_Material_init_Uniform <: TOMLTable
    mean_stress :: Float64
end

mutable struct Input_Material_init_Uniform
    mean_stress :: Float64
end

Base.@kwdef mutable struct TOMLInput_Material_init_K0 <: TOMLTable
    poissons_ratio :: Float64 = NaN
    K0             :: Float64 = isnan(poissons_ratio) ? _undefkey(:K0) : poissons_ratio / (1 - poissons_ratio)
    height_ref     :: Float64
end

mutable struct Input_Material_init_K0
    poissons_ratio :: Float64
    K0             :: Float64
    height_ref     :: Float64
end

#############
# RigidBody #
#############

Base.@kwdef mutable struct TOMLInput_RigidBody <: TOMLTable
    control          :: Bool            = true
    density          :: Float64         = control ? Inf : _undefkey(:density)
    velocity         :: ToVec           = nothing
    angular_velocity :: ToVec           = nothing
    body_force       :: ToVec           = nothing
    model            :: GeometricObject
    function TOMLInput_RigidBody(control, density, velocity, angular_velocity, body_force, model::GeometricObject{dim}) where {dim}
        control                     && (density = Inf)
        isnothing(velocity)         && (velocity = zeros(dim))
        isnothing(angular_velocity) && (angular_velocity = zeros(3))
        isnothing(body_force)       && (body_force = zeros(dim))
        model.m = density * area(model) # TODO: use volume for 3D
        model.v = velocity
        model.ω = angular_velocity
        new(control, density, velocity, angular_velocity, body_force, model)
    end
    function TOMLInput_RigidBody(control, density, velocity, angular_velocity, body_force, model)
        TOMLInput_RigidBody(control, density, velocity, angular_velocity, body_force, convert(GeometricObject, model))
    end
end

mutable struct Input_RigidBody{dim, Model <: GeometricObject{dim}}
    control          :: Bool
    density          :: Float64
    velocity         :: Vec{dim, Float64}
    angular_velocity :: Vec{3, Float64}
    body_force       :: Vec{dim, Float64}
    model            :: Model
end

# Polygon

Base.@kwdef mutable struct TOMLInput_RigidBody_model_Polygon <: TOMLTable
    coordinates :: Vector{Vector{Float64}}
end

function Base.convert(::Type{GeometricObject}, model::TOMLInput_RigidBody_model_Polygon)
    GeometricObject(Polygon(Vec{2}.(model.coordinates)...))
end

# Circle

Base.@kwdef mutable struct TOMLInput_RigidBody_model_Circle <: TOMLTable
    centroid :: Vector{Float64}
    radius   :: Float64
end

function Base.convert(::Type{GeometricObject}, model::TOMLInput_RigidBody_model_Circle)
    GeometricObject(Circle(Vec{2}(model.centroid), model.radius))
end

############
# Advanced #
############

Base.@kwdef mutable struct TOMLInput_Advanced <: TOMLTable
    npoints_in_cell           :: Int     = 2
    contact_threshold_scale   :: Float64 = 1.0
    contact_penalty_parameter :: Float64 = 0.0
    reorder_pointstate        :: Bool    = false
end

mutable struct Input_Advanced
    npoints_in_cell           :: Int
    contact_threshold_scale   :: Float64
    contact_penalty_parameter :: Float64
    reorder_pointstate        :: Bool
end


#############
# TOMLInput #
#############

Base.@kwdef mutable struct TOMLInput <: TOMLTable
    project           :: String                                          = "."
    General           :: TOMLInput_General
    Phase             :: Union{TOMLInput_Phase, Vector{TOMLInput_Phase}}
    BoundaryCondition :: TOMLInput_BoundaryCondition                     = TOMLInput_BoundaryCondition()
    Output            :: TOMLInput_Output
    SoilLayer         :: Vector{TOMLInput_SoilLayer}                     = TOMLInput_SoilLayer[]
    Material          :: Vector{TOMLInput_Material}                      = TOMLInput_Material[]
    RigidBody         :: Vector{TOMLInput_RigidBody}                     = TOMLInput_RigidBody[]
    Advanced          :: TOMLInput_Advanced                              = TOMLInput_Advanced()
    Injection         :: Module                                          = Module()
end

function TOMLInput(dict::Dict{String, Any})::TOMLInput
    construct_input(["TOMLInput"], dict)
end

# helper functions to construct `TOMLInput`
construct_input(current::Vector{String}, value::Any) = value
function construct_input(current::Vector{String}, values::Vector)
    if first(values) isa Dict{String, Any}
        map(value -> construct_input(current, value), values)
    else
        values
    end
end
function construct_input(current::Vector{String}, dict::Dict{String, Any})
    maybe_struct = Symbol(join(current, "_"))
    if isdefined(@__MODULE__, maybe_struct) && isstructtype(eval(maybe_struct))
        T = eval(maybe_struct)
        try
            T(; (Symbol(key) => construct_input([current; key], value) for (key, value) in dict)...)
        catch e
            e isa UndefKeywordError && error(join(current, "."), ": ", e)
            rethrow(e)
        end
    elseif length(dict) == 1 # the case for `model.DruckerPrager`
        key, value = only(dict)
        construct_input([current; key], value)
    else
        error(join(current, "."), " is not defined.")
    end
end

function Base.show(io::IO, input::TOMLTable)
    println(io, typeof(input).name.name, ":")
    len = maximum(length ∘ string, propertynames(input))
    list = map(propertynames(input)) do name
        type = replace(string(fieldtype(typeof(input), name)), "PoingrSimulator." => "", "Poingr." => "")
        string(" ", rpad(name, len), " :: ", type)
    end
    print(io, join(list, '\n'))
end


mutable struct Input{General, Phase, BoundaryCondition, Output, SoilLayer <: Vector, Material <: Vector, RigidBody <: Vector, Advanced, Injection}
    project           :: String
    General           :: General
    Phase             :: Phase
    BoundaryCondition :: BoundaryCondition
    Output            :: Output
    SoilLayer         :: SoilLayer
    Material          :: Material
    RigidBody         :: RigidBody
    Advanced          :: Advanced
    Injection         :: Injection
end

function convert_input(input::TOMLInput, ::Val{:Material})
    materials = getproperty(input, :Material)
    isempty(materials) ? convert_input(input.SoilLayer) : convert_input(materials)
end

# helper functions for convert_input
convert_input(x) = x
convert_input(x::Vector) = map(convert_input, x)
@generated function convert_input(table::TOMLTable)
    names = fieldnames(table)
    exps = [:(convert_input(table, $(Val(name)))) for name in names if !(fieldtype(table, name) <: SkipEntry)]
    T = Symbol(replace(string(table.name.name), "TOML" => ""))
    quote
        $T($(exps...))
    end
end

convert_input(input::TOMLTable, ::Val{name}) where {name} = convert_input(getproperty(input, name))

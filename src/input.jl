using Poingr
using Poingr: Interpolation
using GeometricObjects
using TOML

function parse_input(dict::Dict; project = ".", default_outdir = "output.tmp")
    input = convert_input(TOMLInput(dict))
    input.project = project
    if isempty(input.Output.directory)
        input.Output.directory = default_outdir
    end
    input.Output.directory = joinpath(input.project, input.Output.directory)

    # RigidBody
    for rigidbody in input.RigidBody
        @assert length(rigidbody.Phase) == length(input.Phase)
        @assert length(input.Material) == length(rigidbody.frictions)
        model = rigidbody.model
        for friction in rigidbody.frictions
            len = length(friction)
            if model[] isa Polygon
                @assert len == 1 || len == length(model)
            else
                @assert len == 1
            end
        end
    end

    input.General.type.preprocess_input!(input)
    input
end
function parse_input(str::AbstractString; project = ".", default_outdir = "output.tmp")
    parse_input(TOML.parse(str); project, default_outdir)
end
function parse_inputfile(tomlfile::AbstractString)
    @assert isfile(tomlfile) && endswith(tomlfile, ".toml")
    filename = first(splitext(basename(tomlfile)))
    input = parse_input(read(tomlfile, String); project = dirname(tomlfile), default_outdir = string(filename, ".tmp"))
end

abstract type TOMLInputTable end
abstract type InputTable end

# evaluate expression given as string
struct EvalString{T}
    content::T
end
Base.convert(::Type{EvalString{T}}, x::String) where {T} = EvalString{T}(eval(Meta.parse(x)))
convert_input(x::EvalString) = x.content

# convert Vector to Vec
struct ToVec
    content::Vector{Float64}
end
Base.convert(::Type{ToVec}, x::Vector) = ToVec(x)
Base.isempty(x::ToVec) = isempty(x.content)
convert_input(x::ToVec) = Vec{length(x.content), Float64}(x.content)

# throw undef error if `content` is `missing`
struct MaybeUndef{T}
    content::Union{T, Missing}
end
Base.convert(::Type{MaybeUndef{T}}, x) where {T} = MaybeUndef{T}(x)
convert_input(x::MaybeUndef) = x
Base.isassigned(x::MaybeUndef) = !ismissing(x.content)

undefkeyerror(s) = throw(UndefKeywordError(s))
function Base.getproperty(t::InputTable, name::Symbol)
    v = getfield(t, name)
    if v isa MaybeUndef
        ismissing(v.content) && undefkeyerror(name)
        return v.content
    end
    v
end
function Base.isassigned(t::InputTable, name::Symbol)
    v = try
        getfield(t, name)
    catch
        return false
    end
    ifelse(v isa MaybeUndef, isassigned(v), true)
end

###########
# General #
###########

Base.@kwdef mutable struct TOMLInput_General <: TOMLInputTable
    type              :: EvalString{Module}
    coordinate_system :: String
    domain            :: Vector{Vector{Float64}}
    grid_space        :: Float64
    gravity           :: Float64
    interpolation     :: EvalString{Interpolation} = "LinearWLS(QuadraticBSpline())"
    transfer          :: EvalString{Transfer}      = "Transfer()"
    showprogress      :: Bool                      = true
end

mutable struct Input_General{CoordSystem <: CoordinateSystem, Interp <: Interpolation, Trans <: Transfer} <: InputTable
    type              :: Module
    coordinate_system :: CoordSystem
    domain            :: Vector{Vector{Float64}}
    grid_space        :: Float64
    gravity           :: Float64
    interpolation     :: Interp
    transfer          :: Trans
    showprogress      :: Bool
end

function convert_input(input::TOMLInput_General, ::Val{:coordinate_system})
    c = input.coordinate_system
    c == "planestrain"  && return PlaneStrain()
    c == "axisymmetric" && return Axisymmetric()
    throw(ArgumentError("wrong `coordinate_system`, got \"$c\", use \"planestrain\" or \"axisymmetric\""))
end


#########
# Phase #
#########

Base.@kwdef mutable struct TOMLInput_Phase <: TOMLInputTable
    time_stop     :: Float64
    CFL           :: Float64 = 0.5
    restart       :: String  = ""
    update_motion :: Bool    = true
end

mutable struct Input_Phase <: InputTable
    time_stop     :: Float64
    CFL           :: Float64
    restart       :: String
    update_motion :: Bool
end

#####################
# BoundaryCondition #
#####################

# dirichlet
Base.@kwdef struct TOMLInput_BoundaryCondition_Dirichlet <: TOMLInputTable
    inbounds :: EvalString{Function}
    velocity :: ToVec
    output   :: Bool = true
end
mutable struct Input_BoundaryCondition_Dirichlet{dim} <: InputTable
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

Base.@kwdef mutable struct TOMLInput_BoundaryCondition <: TOMLInputTable
    top    :: Float64 = 0.0
    bottom :: Float64 = 0.0
    left   :: Float64 = 0.0
    right  :: Float64 = 0.0
    Dirichlet :: Vector{TOMLInput_BoundaryCondition_Dirichlet} = []
end

mutable struct Input_BoundaryCondition{Dirichlet} <: InputTable
    top    :: Contact
    bottom :: Contact
    left   :: Contact
    right  :: Contact
    Dirichlet :: Vector{Dirichlet}
end

function convert_input(input::TOMLInput_BoundaryCondition, ::Val{field}) where {field}
    field != :Dirichlet && return Contact(:friction, getproperty(input, field))
    map(convert_input, input.Dirichlet)
end

##########
# Output #
##########

Base.@kwdef mutable struct TOMLInput_Output <: TOMLInputTable
    time_interval  :: Float64
    directory      :: String  = ""
    snapshots      :: Bool    = false
    snapshot_last  :: Bool    = false
    paraview       :: Bool    = true
    paraview_grid  :: Bool    = false
    copy_inputfile :: Bool    = true
    history        :: Bool    = true # only for `PenetrateIntoGround`
    quickview      :: Bool    = false
end

mutable struct Input_Output <: InputTable
    time_interval  :: Float64
    directory      :: String
    snapshots      :: Bool
    snapshot_last  :: Bool
    paraview       :: Bool
    paraview_grid  :: Bool
    copy_inputfile :: Bool
    history        :: Bool
    quickview      :: Bool
end

############
# Material #
############

wrap(x)::Vector{Float64} = x isa Number ? [x] : x

Base.@kwdef mutable struct TOMLInput_Material <: TOMLInputTable
    region  :: EvalString{Function}
    model   :: MaterialModel
    init
end

mutable struct Input_Material{Model <: MaterialModel, Init} <: InputTable
    region  :: Function
    model   :: Model
    init    :: Init
end

#############
# SoilLayer #
#############

Base.@kwdef mutable struct TOMLInput_SoilLayer <: TOMLInputTable
    thickness      :: Float64
    density        :: Float64
    poissons_ratio :: Float64 = NaN
    K0             :: Float64 = isnan(poissons_ratio) ? undefkeyerror(:K0) : poissons_ratio / (1 - poissons_ratio)
    model          :: MaterialModel
end

mutable struct Input_SoilLayer{Model <: MaterialModel} <: InputTable
    thickness      :: Float64
    density        :: Float64
    poissons_ratio :: Float64
    K0             :: Float64
    model          :: Model
end

##############################
# (Material/SoilLayer).model #
##############################

# newtonian fluid

Base.@kwdef mutable struct TOMLInput_Material_model_NewtonianFluid <: TOMLInputTable
    density_ref    :: Float64
    speed_of_sound :: Float64
    viscosity      :: Float64
    second_coefficient_of_viscosity :: Float64 = -2*viscosity/3
end

function Base.convert(::Type{MaterialModel}, model::TOMLInput_Material_model_NewtonianFluid)
    eos = MorrisWaterEOS(;
        c = model.speed_of_sound,
        ρ_ref = model.density_ref,
    )
    NewtonianFluid(
        eos;
        μ = model.viscosity,
        λ = model.second_coefficient_of_viscosity,
    )
end

# Drucker-Prager model

Base.@kwdef mutable struct TOMLInput_Material_model_DruckerPrager <: TOMLInputTable
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
        ϕ = deg2rad(model.friction_angle),
        ψ = deg2rad(model.dilatancy_angle),
        tensioncutoff = model.tension_cutoff,
    )
end

#############################
# (Material/SoilLayer).init #
#############################

Base.@kwdef mutable struct TOMLInput_Material_init_Uniform <: TOMLInputTable
    density     :: MaybeUndef{Float64} = missing # missing is allowed when using NewtonianFluid
    mean_stress :: Float64
end

mutable struct Input_Material_init_Uniform <: InputTable
    density     :: MaybeUndef{Float64}
    mean_stress :: Float64
end

Base.@kwdef mutable struct TOMLInput_Material_init_K0 <: TOMLInputTable
    density        :: Float64
    poissons_ratio :: Float64 = NaN
    K0             :: Float64 = isnan(poissons_ratio) ? undefkeyerror(:K0) : poissons_ratio / (1 - poissons_ratio)
    height_ref     :: Float64
end

mutable struct Input_Material_init_K0 <: InputTable
    density        :: Float64
    poissons_ratio :: Float64
    K0             :: Float64
    height_ref     :: Float64
end

#############
# RigidBody #
#############

Base.@kwdef mutable struct TOMLInput_RigidBody_FrictionWithMaterial <: TOMLInputTable
    coefficient :: Vector{Float64}
    cohesion    :: Vector{Float64} = []
    function TOMLInput_RigidBody_FrictionWithMaterial(coefficient, cohesion)
        new(wrap(coefficient), wrap(cohesion))
    end
end
function convert_input(input::TOMLInput_RigidBody_FrictionWithMaterial)
    μ = input.coefficient
    c = input.cohesion
    c = isempty(c) ? fill(0.0, length(μ)) : c
    [Vec(x) for x in zip(μ, c)]
end

Base.@kwdef mutable struct TOMLInput_RigidBody_Phase <: TOMLInputTable
    control          :: Bool                            = true
    velocity         :: Union{Nothing, Vector{Float64}} = nothing
    angular_velocity :: Union{Nothing, Vector{Float64}} = nothing
    body_force       :: ToVec = []
end
mutable struct Input_RigidBody_Phase{dim} <: InputTable
    control          :: Bool
    velocity         :: Union{Nothing, Vec{dim, Float64}}
    angular_velocity :: Union{Nothing, Vec{3, Float64}}
    body_force       :: Vec{dim, Float64}
end

Base.@kwdef mutable struct TOMLInput_RigidBody <: TOMLInputTable
    Phase                :: Vector{TOMLInput_RigidBody_Phase}
    density              :: Float64 = all(phase->phase.control, Phase) ? Inf : undefkeyerror(:density)
    model                :: GeometricObject
    FrictionWithMaterial :: Vector{TOMLInput_RigidBody_FrictionWithMaterial}
    output               :: Bool = true
    reset_position       :: Bool = true # for PenetrateIntoGround
    # dummies to call convert_input
    control    :: Nothing = nothing
    body_force :: Nothing = nothing
end

mutable struct Input_RigidBody{dim, Model <: GeometricObject{dim}} <: InputTable
    Phase          :: Vector{Input_RigidBody_Phase{dim}}
    density        :: Float64
    model          :: Model
    frictions      :: Vector{Vector{Vec{2, Float64}}}
    output         :: Bool
    reset_position :: Bool
    # to store current phase (TODO: better way)
    control    :: Bool
    body_force :: Vec{dim, Float64}
end

_dimension(::GeometricObject{dim}) where {dim} = dim
function convert_input(input::TOMLInput_RigidBody, ::Val{:Phase})
    dim = _dimension(input.model)
    map(input.Phase) do phase
        Input_RigidBody_Phase{dim}(
            phase.control,
            isnothing(phase.velocity)         ? nothing : convert_input(ToVec(phase.velocity)),
            isnothing(phase.angular_velocity) ? nothing : convert_input(ToVec(phase.angular_velocity)),
            isempty(phase.body_force)         ? zero(Vec{dim}) : convert_input(phase.body_force),
        )
    end
end

convert_input(input::TOMLInput_RigidBody, ::Val{:control}) = true
function convert_input(input::TOMLInput_RigidBody, ::Val{:body_force})
    dim = _dimension(input.model)
    zero(Vec{dim})
end

# Polygon
Base.@kwdef mutable struct TOMLInput_RigidBody_model_Polygon <: TOMLInputTable
    coordinates :: Vector{Vector{Float64}}
end
function Base.convert(::Type{GeometricObject}, model::TOMLInput_RigidBody_model_Polygon)
    GeometricObject(Polygon(Vec{2}.(model.coordinates)...))
end

# Square
Base.@kwdef mutable struct TOMLInput_RigidBody_model_Square <: TOMLInputTable
    centroid :: Vector{Float64}
    radius   :: Float64
    angle    :: Float64         = 0.0
end
function Base.convert(::Type{GeometricObject}, model::TOMLInput_RigidBody_model_Square)
    centroid = model.centroid
    radius = model.radius
    angle = model.angle
    d = radius / √2
    corner1 = centroid + Vec(-d, -d)
    corner2 = centroid + Vec( d, -d)
    corner3 = centroid + Vec( d,  d)
    corner4 = centroid + Vec(-d,  d)
    GeometricObject(rotate(Polygon(corner1, corner2, corner3, corner4), deg2rad(angle)))
end

# Triangle
Base.@kwdef mutable struct TOMLInput_RigidBody_model_Triangle <: TOMLInputTable
    centroid :: Vector{Float64}
    radius   :: Float64
    angle    :: Float64         = 0.0
end
function Base.convert(::Type{GeometricObject}, model::TOMLInput_RigidBody_model_Triangle)
    centroid = model.centroid
    radius = model.radius
    angle = model.angle
    x1 = centroid + radius * Vec(-√3/2, -1/2)
    x2 = centroid + radius * Vec( √3/2, -1/2)
    x3 = centroid + Vec(0, radius)
    GeometricObject(rotate(Polygon(x1, x2, x3), deg2rad(angle)))
end

# Circle
Base.@kwdef mutable struct TOMLInput_RigidBody_model_Circle <: TOMLInputTable
    centroid :: Vector{Float64}
    radius   :: Float64
end
function Base.convert(::Type{GeometricObject}, model::TOMLInput_RigidBody_model_Circle)
    GeometricObject(Circle(Vec{2}(model.centroid), model.radius))
end

############
# Advanced #
############

Base.@kwdef mutable struct TOMLInput_Advanced <: TOMLInputTable
    npoints_in_cell               :: Int     = 2
    contact_threshold_scale       :: Float64 = 1.0
    contact_penalty_parameter     :: Float64 = 0.0
    reorder_pointstate            :: Bool    = false
    dem_contact_penalty_parameter :: Float64 = 0.9
end

mutable struct Input_Advanced <: InputTable
    npoints_in_cell               :: Int
    contact_threshold_scale       :: Float64
    contact_penalty_parameter     :: Float64
    reorder_pointstate            :: Bool
    dem_contact_penalty_parameter :: Float64
end


#############
# TOMLInput #
#############

Base.@kwdef mutable struct TOMLInput <: TOMLInputTable
    project           :: String                      = "."
    General           :: TOMLInput_General
    Phase             :: Vector{TOMLInput_Phase}
    BoundaryCondition :: TOMLInput_BoundaryCondition = TOMLInput_BoundaryCondition()
    Output            :: TOMLInput_Output
    SoilLayer         :: Vector{TOMLInput_SoilLayer} = TOMLInput_SoilLayer[]
    Material          :: Vector{TOMLInput_Material}  = TOMLInput_Material[]
    RigidBody         :: Vector{TOMLInput_RigidBody} = TOMLInput_RigidBody[]
    Advanced          :: TOMLInput_Advanced          = TOMLInput_Advanced()
    Injection         :: Module                      = Module()
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

function Base.show(io::IO, input::TOMLInputTable)
    println(io, typeof(input).name.name, ":")
    len = maximum(length ∘ string, propertynames(input))
    list = map(propertynames(input)) do name
        type = replace(string(fieldtype(typeof(input), name)), "PoingrSimulator." => "", "Poingr." => "")
        string(" ", rpad(name, len), " :: ", type)
    end
    print(io, join(list, '\n'))
end


mutable struct Input{General, Phase, BoundaryCondition, Output, SoilLayer <: Vector, Material <: Vector, RigidBody <: Vector, Advanced, Injection} <: InputTable
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
@generated function convert_input(table::TOMLInputTable)
    names = fieldnames(table)
    exps = [:(convert_input(table, $(Val(name)))) for name in names]
    T = Symbol(replace(string(table.name.name), "TOML" => ""))
    quote
        $T($(exps...))
    end
end

convert_input(input::TOMLInputTable, ::Val{field}) where {field} = convert_input(getproperty(input, field))

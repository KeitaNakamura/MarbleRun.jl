using TOMLX

##########
# Errors #
##########

abstract type InputError <: Exception end

struct UndefKeyError <: InputError
    type::Type
    name::Symbol
end
function Base.showerror(io::IO, e::UndefKeyError)
    print(io, "UndefKeyError: ")
    print(io, "key \"", e.name, "\" must be given for \"", getkeyname(e.type), "\"")
end

struct UnsupportedKeyError <: InputError
    type::Type
    name::Symbol
end
function Base.showerror(io::IO, e::UnsupportedKeyError)
    print(io, "UnsupportedKeyError: ")
    print(io, "got unsupported key \"", e.name, "\" for \"", getkeyname(e.type), "\"")
end

function getkeyname(::Type{T}) where {T}
    name = string(T)
    m = match(r"{.*}", name)
    if m !== nothing
        name = replace(name, m.match=>"")
    end
    replace(name, "MarbleRun."=>"", "Input_"=>"", "_"=>".")
end

#############
# readinput #
#############

function readinput(dict::AbstractDict; project = ".", default_outdir = "output.tmp")
    # read input
    if dict[:General][:type] == FreeRun
        InputType = Input{2, Input_Material}
    elseif dict[:General][:type] == GroundPenetration
        InputType = Input{2, Input_SoilLayer}
    else
        error()
    end
    input = try
        TOMLX.from_dict(InputType, dict)
    catch e
        if e isa TOMLX.FieldError
            e isa TOMLX.UndefFieldError       && throw(UndefKeyError(e.type, e.name))
            e isa TOMLX.UnsupportedFieldError && throw(UnsupportedKeyError(e.type, e.name))
            rethrow()
        else
            rethrow()
        end
    end

    # check version
    input_version = VersionNumber(input.MarbleRun)
    if Base.thisminor(input_version) != Base.thisminor(PKG_VERSION)
        @warn """
        Current MarbleRun version $PKG_VERSION is different from version $input_version in input file.
        Unexpected behavior may occur.
        """
    end

    # project/output directory
    input.project = project
    if isempty(input.Output.directory)
        input.Output.directory = default_outdir
    end
    input.Output.directory = joinpath(input.project, input.Output.directory)

    # quickview
    if input.Output.quickview
        if input.General.showprogress == false
            @info "force `showprogress=true` by `quickview=true`"
            input.General.showprogress = true
        end
    end

    # RigidBody
    for rigidbody in input.RigidBody
        @assert length(rigidbody.Phase) == length(input.Phase)
        @assert length(input.Material) == length(rigidbody.FrictionWithMaterial)
        model = rigidbody.model
        for coef in rigidbody.FrictionWithMaterial
            len = length(coef)
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
function readinput(str::AbstractString; project = ".", default_outdir = "output.tmp")
    readinput(TOMLX.parse(@__MODULE__, str); project, default_outdir)
end
function readinputfile(tomlfile::AbstractString)
    @assert isfile(tomlfile) && endswith(tomlfile, ".toml")
    filename = first(splitext(basename(tomlfile)))
    input = readinput(read(tomlfile, String); project = dirname(tomlfile), default_outdir = string(filename, ".tmp"))
end

undefkeyerror(::Type{T}, name::Symbol) where {T} = throw(UndefKeyError(T, name))

###########
# General #
###########

Base.@kwdef mutable struct Input_General
    type              :: Module
    coordinate_system :: CoordinateSystem
    domain            :: Vector{Vector{Float64}}
    grid_space        :: Float64
    gravity           :: Float64
    interpolation     :: Interpolation = LinearWLS(QuadraticBSpline())
    transfer          :: Transfer      = Transfer(interpolation)
    v_p_formulation   :: Bool          = false
    showprogress      :: Bool          = true
end

#########
# Phase #
#########

Base.@kwdef struct Input_Phase
    time_stop     :: Float64
    CFL           :: Float64
    restart       :: String  = ""
    update_motion :: Bool    = true
end

#####################
# BoundaryCondition #
#####################

# dirichlet
Base.@kwdef mutable struct Input_BoundaryCondition_Dirichlet{dim}
    region         :: Function
    velocity       :: Vec{dim, Float64}
    output         :: Bool = true
    # private
    displacement   :: Float64 = 0.0
    reaction_force :: Float64 = 0.0
    node_indices   :: Vector{CartesianIndex{dim}} = CartesianIndex{dim}[]
end

Base.@kwdef mutable struct Input_BoundaryCondition{dim}
    top    :: Float64 = 0.0
    bottom :: Float64 = 0.0
    left   :: Float64 = 0.0
    right  :: Float64 = 0.0
    sides  :: Vector{Pair{String, CoulombFriction}} = [ # private
        "-x" => CoulombFriction(; μ=left),
        "+x" => CoulombFriction(; μ=right),
        "-y" => CoulombFriction(; μ=bottom),
        "+y" => CoulombFriction(; μ=top),
    ]
    Dirichlet :: Vector{Input_BoundaryCondition_Dirichlet{dim}} = []
end

##########
# Output #
##########

Base.@kwdef mutable struct Input_Output
    time_interval  :: Float64
    directory      :: String  = ""
    snapshots      :: Bool    = false
    snapshot_first :: Bool    = false
    snapshot_last  :: Bool    = false
    copy_inputfile :: Bool    = true
    history        :: Bool    = true  # only for `GroundPenetration`
    quickview      :: Bool    = false
end

############
# Paraview #
############

Base.@kwdef struct Input_Paraview_PointState
    velocity          :: Bool = false
    displacement      :: Bool = false
    mean_stress       :: Bool = false
    pressure          :: Bool = false
    deviatoric_stress :: Bool = false
    volumetric_strain :: Bool = false
    deviatoric_strain :: Bool = false
    stress            :: Bool = false
    strain            :: Bool = false
    vorticity         :: Bool = false
    density           :: Bool = false
    material_index    :: Bool = false
end

Base.@kwdef struct Input_Paraview_GridState
    velocity         :: Bool = false
    contact_force    :: Bool = false
    contact_distance :: Bool = false
end

Base.@kwdef mutable struct Input_Paraview
    output     :: Bool                                     = true
    PointState :: Input_Paraview_PointState                = Input_Paraview_PointState(; velocity = true)
    GridState  :: Union{Nothing, Input_Paraview_GridState} = nothing
end

############
# Material #
############

# Init
abstract type Init end

Base.@kwdef struct InitUniform <: Init
    density     :: Union{Nothing, Float64} = nothing # nothing is allowed when using NewtonianFluid
    mean_stress :: Float64
end

Base.@kwdef struct InitK0 <: Init
    density        :: Float64
    poissons_ratio :: Float64 = NaN
    K0             :: Float64 = isnan(poissons_ratio) ? undefkeyerror(InitK0, :K0) : poissons_ratio / (1 - poissons_ratio)
    height_ref     :: Float64
end

wrap(x)::Vector{Float64} = x isa Number ? [x] : x

Base.@kwdef struct Input_Material
    region :: Union{Function, String}
    model  :: MaterialModel
    init   :: Init
end

#############
# SoilLayer #
#############

Base.@kwdef struct Input_SoilLayer
    thickness      :: Float64
    density        :: Float64
    poissons_ratio :: Float64 = NaN
    K0             :: Float64 = isnan(poissons_ratio) ? (undefkeyerror(Input_SoilLayer, :K0)) : poissons_ratio / (1 - poissons_ratio)
    model          :: MaterialModel
end

###################
# Material models #
###################

Base.@kwdef struct MaterialNewtonianFluid
    density_ref                     :: Float64
    speed_of_sound                  :: Float64
    dynamic_viscosity               :: Float64
    second_coefficient_of_viscosity :: Float64 = -2*dynamic_viscosity/3
end

function Base.convert(::Type{MaterialModel}, model::MaterialNewtonianFluid)
    eos = MorrisWaterEOS(;
        c = model.speed_of_sound,
        ρ_ref = model.density_ref,
    )
    NewtonianFluid(
        eos;
        μ = model.dynamic_viscosity,
        λ = model.second_coefficient_of_viscosity,
    )
end

Base.@kwdef struct MaterialDruckerPrager
    mohr_coulomb_type :: String
    youngs_modulus    :: Float64
    poissons_ratio    :: Float64
    cohesion          :: Float64
    friction_angle    :: Float64
    dilatancy_angle   :: Float64
    tension_cutoff    :: Union{Float64, Bool} = false
end

function Base.convert(::Type{MaterialModel}, model::MaterialDruckerPrager)
    DruckerPrager(
        LinearElastic(; E = model.youngs_modulus, ν = model.poissons_ratio),
        model.mohr_coulomb_type;
        c = model.cohesion,
        ϕ = deg2rad(model.friction_angle),
        ψ = deg2rad(model.dilatancy_angle),
        tensioncutoff = model.tension_cutoff,
    )
end

#############
# RigidBody #
#############

struct Input_RigidBody_FrictionWithMaterial <: AbstractVector{Vec{2, Float64}}
    coefs::Vector{Vec{2, Float64}}
end
Base.size(x::Input_RigidBody_FrictionWithMaterial) = size(x.coefs)
Base.getindex(x::Input_RigidBody_FrictionWithMaterial, i::Int) = (@_propagate_inbounds_meta; x.coefs[i])

function Input_RigidBody_FrictionWithMaterial(; coefficient, cohesion=[])
    μ = wrap(coefficient)
    c = wrap(cohesion)
    c = isempty(c) ? fill(0.0, length(μ)) : c
    Input_RigidBody_FrictionWithMaterial([Vec(x) for x in zip(μ, c)])
end

Base.@kwdef struct Input_RigidBody_Phase{dim}
    control          :: Bool                              = true
    velocity         :: Union{Nothing, Vec{dim, Float64}} = nothing
    angular_velocity :: Union{Nothing, Vec{3, Float64}}   = nothing
    body_force       :: Vec{dim, Float64}                 = zero(Vec{dim})
end

Base.@kwdef mutable struct Input_RigidBody{dim}
    model                :: GeometricObject{dim, Float64}
    Phase                :: Vector{Input_RigidBody_Phase{dim}}
    FrictionWithMaterial :: Vector{Input_RigidBody_FrictionWithMaterial}
    density              :: Float64                                     = all(phase->phase.control, Phase) ? Inf : undefkeyerror(Input_RigidBody, :density)
    inverse              :: Bool                                        = false
    output               :: Bool                                        = true
    reset_position       :: Bool                                        = true # for GroundPenetration
    # to store current phase (TODO: better way)
    control    :: Bool                                                  = true
    body_force :: Vec{dim, Float64}                                     = zero(Vec{dim})
end

# Polygon
Base.@kwdef struct RigidPolygon
    coordinates :: Vector{Vector{Float64}}
end
function Base.convert(::Type{<: GeometricObject{2}}, model::RigidPolygon)
    GeometricObject(Polygon(Vec{2}.(model.coordinates)...))
end

# Square
Base.@kwdef struct RigidSquare
    centroid :: Vector{Float64}
    radius   :: Float64
    angle    :: Float64 = 0.0
end
function Base.convert(::Type{<: GeometricObject{2}}, model::RigidSquare)
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
Base.@kwdef struct RigidTriangle
    centroid :: Vector{Float64}
    radius   :: Float64
    angle    :: Float64 = 0.0
end
function Base.convert(::Type{<: GeometricObject{2}}, model::RigidTriangle)
    centroid = model.centroid
    radius = model.radius
    angle = model.angle
    x1 = centroid + radius * Vec(-√3/2, -1/2)
    x2 = centroid + radius * Vec( √3/2, -1/2)
    x3 = centroid + Vec(0, radius)
    GeometricObject(rotate(Polygon(x1, x2, x3), deg2rad(angle)))
end

# Circle
Base.@kwdef struct RigidCircle
    centroid :: Vector{Float64}
    radius   :: Float64
end
function Base.convert(::Type{<: GeometricObject{2}}, model::RigidCircle)
    GeometricObject(Circle(Vec{2}(model.centroid), model.radius))
end

############
# Advanced #
############

Base.@kwdef struct Input_Advanced
    npoints_in_cell                           :: Int             = 2
    random_points_generation                  :: Bool            = false
    contact_threshold_scale                   :: Float64         = 1.0
    contact_threshold_scale_for_initial_state :: Float64         = contact_threshold_scale
    contact_penalty_parameter                 :: Float64         = 0.0
    reorder_pointstate                        :: Bool            = false
    dem_contact_penalty_parameter             :: Float64         = 0.9
    dem_contact_for_DKT                       :: Vector{Float64} = Float64[] # [k,ϵ] (experimental)
end

#########
# Input #
#########

Base.@kwdef mutable struct Input{dim, Mat}
    MarbleRun         :: String
    project           :: String                       = "."
    General           :: Input_General
    Phase             :: Vector{Input_Phase}
    BoundaryCondition :: Input_BoundaryCondition{dim} = Input_BoundaryCondition{dim}()
    Output            :: Input_Output
    Paraview          :: Input_Paraview
    SoilLayer         :: Vector{Input_SoilLayer}      = Input_SoilLayer[]
    Material          :: Vector{Mat}                  = SoilLayer
    RigidBody         :: Vector{Input_RigidBody{dim}} = Input_RigidBody{dim}[]
    Advanced          :: Input_Advanced               = Input_Advanced()
    Injection         :: Module                       = Module()
end

module MarbleRun

using Reexport
@reexport using Marble
using MaterialModels
using GeometricObjects

using Serialization
using UnicodePlots
using TOMLX

using Base: @_propagate_inbounds_meta, @_inline_meta

const PKG_VERSION = let
    project = TOMLX.parsefile(joinpath(pkgdir(@__MODULE__), "Project.toml"))
    VersionNumber(project["version"])
end

include("input.jl")
include("methods.jl")
include("dem.jl")
include("GroundPenetration.jl")
include("FreeRun.jl")

# for PackageCompiler
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

##################
# main functions #
##################

# helpers
commas(num::Integer) = replace(string(num), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
function replace_version(toml::AbstractString, version::VersionNumber)
    input_version = TOMLX.parse(@__MODULE__, toml)["MarbleRun"]
    lines = split(toml, '\n')
    n = findfirst(lines) do line
        startswith(replace(line, " " => ""), "MarbleRun=\"$input_version\"")
    end
    lines[n] = "MarbleRun = \"$version\" # \"$input_version\" is replaced by MarbleRun"
    join(lines, '\n')
end

function main(toml::AbstractString; kwargs...)
    istomlfile = try
        isfile(toml) && endswith(toml, ".toml")
    catch e # handle error "name too long" in `isfile`
        false
    end
    if istomlfile
        main_tomlfile(toml; kwargs...)
    else
        main_tomlstring(toml; kwargs...)
    end
end

# from toml file
function main_tomlfile(tomlfile::AbstractString)
    @assert isfile(tomlfile) && endswith(tomlfile, ".toml")
    project = dirname(tomlfile)
    injection_file = joinpath(project, "Injection.jl")
    filename = first(splitext(basename(tomlfile)))
    main_tomlstring(
        read(tomlfile, String);
        injection = isfile(injection_file) ? include(injection_file) : Module(),
        project,
        default_outdir = string(filename, ".tmp"),
    )
end

# from toml string
function main_tomlstring(toml::AbstractString; injection::Module = Module(), project::AbstractString = ".", default_outdir::AbstractString = "output.tmp")
    input = readinput(toml; project, default_outdir)
    input.Injection = injection
    mkpath(input.Output.directory)
    if input.Output.copy_inputfile
        toml = replace_version(toml, PKG_VERSION)
        write(joinpath(input.Output.directory, "input.toml"), toml)
    end
    main(input, input.Phase)
end

# main phases
function main(input::Input, phases::Vector{Input_Phase})
    outdir = input.Output.directory
    t, grid, gridstate, pointstate, rigidbodies, data... = initialize(input, first(phases))
    for i in eachindex(phases)
        input.Output.directory = joinpath(outdir, string(i))
        reinitialize!(rigidbodies, input.RigidBody, i)
        t = main(input, phases[i], (t, grid, gridstate, pointstate, rigidbodies, data...))
    end
end

# main for each phase
function main(input::Input, phase::Input_Phase, (t, grid, gridstate, pointstate, data...) = initialize(input, phase))
    println("Points: ", commas(length(pointstate)))
    input.General.mod.main(input, phase, t, grid, gridstate, pointstate, data...)
end

############################
# initialize/reinitialize! #
############################

function initialize(input::Input, phase::Input_Phase)
    if isempty(phase.restart)
        input.General.mod.initialize(input)
    else
        deserialize(joinpath(input.project, phase.restart))
    end
end

function reinitialize!(rigidbody::GeometricObject, input::Input_RigidBody, phase_index::Int)
    phase = input.Phase[phase_index]
    if phase.control
        rigidbody.m = Inf
    else
        rigidbody.m = input.density * area(geometry(rigidbody)) # TODO: use volume for 3D
    end
    if phase.velocity !== nothing
        rigidbody.v = phase.velocity
    end
    if phase.angular_velocity !== nothing
        rigidbody.ω = phase.angular_velocity
    end
    input.control = phase.control
    input.body_force = phase.body_force
end
function reinitialize!(rigidbodies::Vector{<: GeometricObject}, inputs::Vector{<: Input_RigidBody}, phase_index::Int)
    for (rigidbody, input) in zip(rigidbodies, inputs)
        reinitialize!(rigidbody, input, phase_index)
    end
end
function reinitialize!(rigidbody::GeometricObject, inputs::Vector{<: Input_RigidBody}, phase_index::Int)
    reinitialize!(rigidbody, only(inputs), phase_index)
end
reinitialize!(rigidbodies::Vector, inputs::Vector, phase_index::Int) = @assert isempty(rigidbodies) && isempty(inputs)

end # module

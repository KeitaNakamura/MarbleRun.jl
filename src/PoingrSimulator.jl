module PoingrSimulator

using Poingr
using GeometricObjects

using TOML
using Serialization

using Base: @_propagate_inbounds_meta, @_inline_meta

include("input.jl")
include("methods.jl")
include("dem.jl")
include("PenetrateIntoGround.jl")
include("FreeRun.jl")


function main(tomlfile::AbstractString)
    @assert isfile(tomlfile) && endswith(tomlfile, ".toml")
    project = dirname(tomlfile)
    injection_file = joinpath(project, "injection.jl")
    filename = first(splitext(basename(tomlfile)))
    main(
        read(tomlfile, String),
        isfile(injection_file) ? include(injection_file) : Module();
        project,
        default_outdir = string(filename, ".tmp"),
    )
end

function main(inputtoml::AbstractString, Injection::Module; project::AbstractString = ".", default_outdir::AbstractString = "output.tmp")
    input = parse_input(inputtoml; project, default_outdir)
    input.Injection = Injection
    mkpath(input.Output.directory)
    if input.Output.copy_inputfile
        write(joinpath(input.Output.directory, "input.toml"), inputtoml)
    end
    main(input, input.Phase)
end

function main(input::Input, phase::Input_Phase, (t, grid, pointstate, data...) = initialize(input, phase))
    println("Points: ", length(pointstate))
    @eval $input.General.type.main($input, $phase, $t, $grid, $pointstate, $data...)
end

function main(input::Input, phases::Vector{Input_Phase})
    outdir = input.Output.directory
    t, data... = initialize(input, first(phases))
    for i in eachindex(phases)
        input.Output.directory = joinpath(outdir, string(i))
        t = main(input, phases[i], (t, data...))
    end
end

function initialize(input::Input, phase::Input_Phase)
    if isempty(phase.restart)
        @eval $input.General.type.initialize($input)
    else
        deserialize(joinpath(input.project, phase.restart))
    end
end

end # module

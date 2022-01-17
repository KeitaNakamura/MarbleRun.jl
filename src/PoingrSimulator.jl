module PoingrSimulator

using Poingr
using GeometricObjects
using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

function main(inputtoml_file::AbstractString)
    proj_dir = dirname(inputtoml_file)
    injection_file = joinpath(proj_dir, "injection.jl")
    main(inputtoml_file, isfile(injection_file) ? include(injection_file) : Module())
end

function main(inputtoml_file::AbstractString, Injection::Module)
    proj_dir = dirname(inputtoml_file)
    inputtoml = read(inputtoml_file, String)
    main(proj_dir, inputtoml, Injection)
end

function main(proj_dir::AbstractString, inputtoml::AbstractString, Injection::Module)
    dict = TOML.parse(inputtoml)
    dict["General"]["project_directory"] = proj_dir
    dict["Output"]["directory"] = joinpath(proj_dir, dict["Output"]["directory"])
    dict["Injection"] = Injection
    INPUT = parse_input(dict)

    mkpath(INPUT.Output.directory)
    if INPUT.Output.copy_inputfile
        write(joinpath(INPUT.Output.directory, "input.toml"), inputtoml)
    end

    # use eval for error related with `Injection.main_output`: "method too new to be called from this world context."
    # don't know the mechanism
    @eval $INPUT.General.simulation.main($INPUT)
end

include("input.jl")
include("methods.jl")
include("PenetrateIntoGround.jl")
include("FreeRun.jl")

end # module

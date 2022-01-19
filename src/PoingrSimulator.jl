module PoingrSimulator

using Poingr
using GeometricObjects
using TOML
using Serialization

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

function build_INPUT(proj_dir::AbstractString, inputtoml::AbstractString, Injection::Module)
    dict = TOML.parse(inputtoml)
    dict["General"]["project_directory"] = proj_dir
    dict["Injection"] = Injection
    if haskey(dict["General"], "restart")
        dir_ext = splitext(dict["Output"]["directory"]) # handle .tmp folder
        output_dir = string(dir_ext[1], "_restarted_from_", dict["General"]["restart"], dir_ext[2])
        dict["Output"]["original_directory"] = joinpath(proj_dir, dict["Output"]["directory"])
        dict["Output"]["directory"] = joinpath(proj_dir, output_dir)
    else
        dict["Output"]["directory"] = joinpath(proj_dir, dict["Output"]["directory"])
    end
    INPUT = parse_input(dict)
    INPUT
end

function build_INPUT(inputtoml_file::AbstractString)
    @assert isfile(inputtoml_file)
    proj_dir = dirname(inputtoml_file)
    inputtoml = read(inputtoml_file, String)
    injection_file = joinpath(proj_dir, "injection.jl")
    build_INPUT(proj_dir, inputtoml, isfile(injection_file) ? include(injection_file) : Module())
end

function main(proj_dir::AbstractString, inputtoml::AbstractString, Injection::Module)
    INPUT = build_INPUT(proj_dir, inputtoml, Injection)

    mkpath(INPUT.Output.directory)
    if INPUT.Output.copy_inputfile
        write(joinpath(INPUT.Output.directory, "input.toml"), inputtoml)
    end

    # use eval for error related with `Injection.main_output`: "method too new to be called from this world context."
    # don't know the mechanism
    if haskey(INPUT.General, :restart)
        index = string(INPUT.General.restart)
        data = deserialize(joinpath(INPUT.Output.original_directory, "snapshots", "snapshot$index"))
        @eval $INPUT.General.type.main($INPUT, $data...)
    else
        @eval $INPUT.General.type.main($INPUT, $INPUT.General.type.initialize($INPUT)...)
    end
end

include("input.jl")
include("methods.jl")
include("PenetrateIntoGround.jl")
include("FreeRun.jl")

end # module

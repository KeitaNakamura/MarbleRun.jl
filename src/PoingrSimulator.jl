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
    INPUT = parse_input(dict)

    # create output directory
    output_dir = joinpath(proj_dir, INPUT.Output.folder_name)
    mkpath(output_dir)

    # copy input toml file
    if INPUT.Output.copy_inputfile
        write(joinpath(output_dir, "input.toml"), inputtoml)
    end

    simulation = Symbol(INPUT.General.simulation)
    @eval $simulation.main($proj_dir, $INPUT, $Injection)
end

include("input.jl")
include("methods.jl")
include("PenetrateIntoGround.jl")
include("FreeRun.jl")

end # module

module PoingrSimulator

using Poingr
using GeometricObjects
using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

function main(inputtoml_file::AbstractString)
    proj_dir = splitdir(inputtoml_file)[1]
    inputtoml = read(inputtoml_file, String)
    main(proj_dir, inputtoml)
end

function main(proj_dir::AbstractString, inputtoml::AbstractString)
    INPUT = parseinput(TOML.parse(inputtoml))

    # create output directory
    output_dir = joinpath(proj_dir, INPUT.Output.folder_name)
    mkpath(output_dir)

    # copy input toml file
    if INPUT.Output.copy_inputfile
        write(joinpath(output_dir, "input.toml"), inputtoml)
    end

    # custom output
    custom_output_file = joinpath(proj_dir, INPUT.Output.custom_output)
    if isfile(custom_output_file)
        CustomOutput = include(custom_output_file)
        @assert CustomOutput isa Module && isdefined(CustomOutput, :main)
        INPUT = (; CustomOutput, INPUT...)
    end

    simulation = Symbol(INPUT.General.simulation)
    @eval $simulation.main($proj_dir, $INPUT)
end

include("utils.jl")
include("PenetrateIntoGround.jl")

end # module

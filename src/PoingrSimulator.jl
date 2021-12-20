module PoingrSimulator

using Poingr
using GeometricObjects
using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

function main(inputtoml_file::AbstractString, Injection::Module = Module())
    proj_dir = splitdir(inputtoml_file)[1]
    inputtoml = read(inputtoml_file, String)
    main(proj_dir, inputtoml, Injection)
end

function main(proj_dir::AbstractString, inputtoml::AbstractString, Injection::Module)
    INPUT = parseinput(TOML.parse(inputtoml))

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

include("utils.jl")
include("PenetrateIntoGround.jl")

end # module

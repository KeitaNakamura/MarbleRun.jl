module PoingrSimulator

using Poingr
using GeometricObjects
using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

function main(inputtoml::AbstractString)
    proj_dir = splitdir(inputtoml)[1]
    INPUT = parseinput(inputtoml)

    # create output directory
    output_dir = joinpath(proj_dir, INPUT.Output.folder_name)
    mkpath(output_dir)

    # copy input toml file
    if INPUT.Output.copy_inputfile
        cp(inputtoml, joinpath(output_dir, "input.toml"); force = true)
    end

    simulation = Symbol(INPUT.General.simulation)
    @eval $simulation.main($proj_dir, $INPUT)
end

include("utils.jl")
include("PenetrateIntoGround.jl")

end # module

module PoingrSimulator

using Poingr
using GeometricObjects
using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

function preprocess_input!(dict::Dict)
    if haskey(dict["General"], "coordinate_system")
        coordinate_system = dict["General"]["coordinate_system"]
        if coordinate_system == "plane_strain"
            dict["General"]["coordinate_system"] = PlaneStrain()
        elseif coordinate_system == "axisymmetric"
            dict["General"]["coordinate_system"] = Axisymmetric()
        else
            throw(ArgumentError("wrong `coordinate_system`, got \"$coordinate_system\", use \"plane_strain\" or \"axisymmetric\""))
        end
    end
    if haskey(dict, "Material")
        for mat in dict["Material"]
            if haskey(mat, "region")
                mat["region"] = eval(Meta.parse(mat["region"])) # should be anonymous function
            end
            if haskey(mat, "type")
                mat["type"] = eval(Meta.parse(mat["type"]))
            end
        end
    end
end

function main(inputtoml_file::AbstractString, Injection::Module = Module())
    proj_dir = splitdir(inputtoml_file)[1]
    inputtoml = read(inputtoml_file, String)
    main(proj_dir, inputtoml, Injection)
end

function main(proj_dir::AbstractString, inputtoml::AbstractString, Injection::Module)
    dict = TOML.parse(inputtoml)
    preprocess_input!(dict)
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

include("utils.jl")
include("transfer.jl")
include("PenetrateIntoGround.jl")
include("FreeRun.jl")

end # module

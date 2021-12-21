using PoingrSimulator
using Test

using TOML
using CSV

# used in injection.jl
using DelimitedFiles

const fix_results = false

function check_results(inputtoml::String)
    @assert endswith(inputtoml, ".toml")
    testname = first(splitext(basename(inputtoml)))
    @testset "$(joinpath(basename(dirname(inputtoml)), basename(inputtoml)))" begin
        injection_file = joinpath(dirname(inputtoml), "injection.jl")
        if isfile(injection_file)
            PoingrSimulator.main(inputtoml, include(injection_file))
        else
            PoingrSimulator.main(inputtoml)
        end
        println()
        output_dir = TOML.parsefile(inputtoml)["Output"]["folder_name"]

        if fix_results
            mv(joinpath(dirname(inputtoml), output_dir, "history.csv"),
               joinpath(dirname(inputtoml), "output", "$testname.csv"); force = true)
        else
            # check results
            output = CSV.File(joinpath(dirname(inputtoml), "output", "$testname.csv")) # expected output
            history = CSV.File(joinpath(dirname(inputtoml), output_dir, "history.csv"))
            for name in propertynames(output)
                output_col = output[name]
                history_col = history[name]
                @test output_col ≈ history_col  rtol = 1e-3
            end
        end
    end
end

@testset "PenetrateIntoGround" begin
    for (root, dirs, files) in walkdir("PenetrateIntoGround")
        for file in files
            path = joinpath(root, file)
            basename(dirname(path)) != "output.tmp" && endswith(path, ".toml") && check_results(path)
        end
    end
end

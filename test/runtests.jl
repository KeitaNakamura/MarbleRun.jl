using PoingrSimulator
using Test

using TOML
using CSV

function check_results(inputtoml::String)
    @assert endswith(inputtoml, ".toml")
    testname = first(splitext(basename(inputtoml)))
    @testset "$testname" begin
        PenetrateIntoGround.main(inputtoml); println()
        output_dir = TOML.parsefile(inputtoml)["Output"]["folder_name"]

        # check results
        output = CSV.File(joinpath(dirname(inputtoml), "output", "$testname.csv")) # expected output
        history = CSV.File(joinpath(dirname(inputtoml), output_dir, "history.csv"))
        for name in propertynames(output)
            output_col = output[name]
            history_col = history[name]
            for i in 1:length(output_col)
                val = output_col[i]
                @test 0.98*val ≤ history_col[i] ≤ 1.02*val # ±2%
            end
        end
    end
end

@testset "PenetrateIntoGround" begin
    for (root, dirs, files) in walkdir("PenetrateIntoGround")
        for file in files
            path = joinpath(root, file)
            endswith(path, ".toml") && check_results(path)
        end
    end
end

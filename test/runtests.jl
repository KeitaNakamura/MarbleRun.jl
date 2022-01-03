using PoingrSimulator
using Test

using TOML
using CSV

using LinearAlgebra # norm
using ReadVTK
using NaturalSort

# used in injection.jl
using DelimitedFiles

const fix_results = false

function check_results(inputtoml::String; check_history = false)
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

        INPUT = PoingrSimulator.parseinput(TOML.parsefile(inputtoml))
        proj_dir = dirname(inputtoml)
        output_dir = INPUT.Output.folder_name

        vtk_file = joinpath(
            "paraview",
            sort(
                filter(
                    file -> endswith(file, ".vtu"),
                    only(walkdir(joinpath(proj_dir, output_dir, "paraview")))[3]
                ),
                lt = natural
            )[end],
        )

        if fix_results
            cp(joinpath(proj_dir, output_dir, vtk_file),
               joinpath(proj_dir, "output", "$testname.vtu"); force = true)
        else
            # check results
            expected = VTKFile(joinpath(proj_dir, "output", "$testname.vtu")) # expected output
            result = VTKFile(joinpath(proj_dir, output_dir, vtk_file))
            expected_points = get_points(expected)
            result_points = get_points(result)
            @assert size(expected_points) == size(result_points)
            @test all(eachindex(expected_points)) do i
                norm(expected_points[i] - result_points[i]) < 0.05*INPUT.General.grid_space
            end
        end

        if check_history
            if fix_results
                cp(joinpath(proj_dir, output_dir, "history.csv"),
                   joinpath(proj_dir, "output", "$testname.csv"); force = true)
            else
                # check results
                output = CSV.File(joinpath(proj_dir, "output", "$testname.csv")) # expected output
                history = CSV.File(joinpath(proj_dir, output_dir, "history.csv"))
                for name in propertynames(output)
                    output_col = output[name]
                    history_col = history[name]
                    @test output_col ≈ history_col  rtol = 1e-3
                end
            end
        end
    end
end

@testset "PenetrateIntoGround" begin
    for (root, dirs, files) in walkdir("PenetrateIntoGround")
        for file in files
            path = joinpath(root, file)
            basename(dirname(path)) != "output.tmp" && endswith(path, ".toml") &&
                check_results(path; check_history = true)
        end
    end
end

using PoingrSimulator
using Test

using LinearAlgebra # norm
using CSV           # history file
using Serialization # snapshots

# vtk
using ReadVTK
using NaturalSort

const fix_results = false

function check_results(tomlfile::String)
    @assert endswith(tomlfile, ".toml")
    testcase = joinpath(basename(dirname(tomlfile)), basename(tomlfile))
    println("\n>> ", testcase)
    @testset "$testcase" begin
        @time PoingrSimulator.main(tomlfile)

        input = PoingrSimulator.parse_inputfile(tomlfile)
        for phase_index in 1:length(input.Phase)
            proj_dir = input.project
            output_dir = joinpath(input.Output.directory, string(phase_index))

            restart_case = !isempty(input.Phase[phase_index].restart)

            # for restart
            testname = first(splitext(basename(tomlfile)))
            if restart_case
                testname = replace(testname, "_restart" => "")
            end
            testname *= "_$phase_index"

            # vtk files
            vtk_file = joinpath(
                "paraview",
                sort(
                    filter(
                        file -> endswith(file, "_1.vtu"),
                        only(walkdir(joinpath(output_dir, "paraview")))[3]
                    ),
                    lt = natural
                )[end],
            )

            if fix_results
                cp(joinpath(output_dir, vtk_file),
                   joinpath(proj_dir, "output", "$testname.vtu"); force = true)
            elseif !restart_case
                # check results
                expected = VTKFile(joinpath(proj_dir, "output", "$testname.vtu")) # expected output
                result = VTKFile(joinpath(output_dir, vtk_file))
                expected_points = get_points(expected)
                result_points = get_points(result)
                @test size(expected_points) == size(result_points)
                @test all(eachindex(expected_points)) do i
                    norm(expected_points[i] - result_points[i]) < 0.05*input.General.grid_space
                end
            end

            # snapshots file
            if !fix_results && input.Output.snapshots == true
                root, _, files = only(walkdir(joinpath(output_dir, "snapshots")))
                count = 0
                for file in sort(files, lt = natural)
                    @test file == "snapshot$count"
                    @test deserialize(joinpath(root, file)) isa NamedTuple
                    count += 1
                end
            end

            # test history.csv if exists
            if isfile(joinpath(output_dir, "history.csv")) && !restart_case
                expected = joinpath(output_dir, "history.csv")
                src = joinpath(proj_dir, "output", "$testname.csv")
                fix_results ? fix_history(expected, src) : check_history(expected, src)
            end

            # dirichlet
            if any(d -> d.output, input.BoundaryCondition.Dirichlet)
                Dirichlet = input.BoundaryCondition.Dirichlet
                for i in eachindex(Dirichlet)
                    dirichlet = Dirichlet[i]
                    if dirichlet.output
                        expected = joinpath(output_dir, "dirichlet", "$i", "history.csv")
                        src = joinpath(proj_dir, "output", "$(testname)_diriclet_$i.csv")
                        fix_results ? fix_history(expected, src) : check_history(expected, src)
                    end
                end
            end

            # rigidbodies
            if any(d -> d.output, input.RigidBody) && startswith(tomlfile, "FreeRun")
                RigidBody = input.RigidBody
                for i in eachindex(RigidBody)
                    rigidbody = RigidBody[i]
                    if rigidbody.output
                        expected = joinpath(output_dir, "rigidbodies", "$i", "history.csv")
                        src = joinpath(proj_dir, "output", "$(testname)_rigidbody_$i.csv")
                        fix_results ? fix_history(expected, src) : check_history(expected, src)
                    end
                end
            end
        end
    end
end

function fix_history(dest, src)
    cp(dest, src ; force = true)
end
function check_history(expected, src)
    # check results
    output = CSV.File(expected) # expected output
    history = CSV.File(src)
    for name in propertynames(output)
        output_col = output[name]
        history_col = history[name]
        @test output_col ≈ history_col atol=1e-8 rtol=1e-3
    end
end

@testset "$module_name" for module_name in ("PenetrateIntoGround", "FreeRun",)
    # clean up  first
    for (root, dirs, files) in collect(walkdir(module_name))
        for dir in dirs
            endswith(dir, ".tmp") && rm(joinpath(root, dir); recursive = true)
        end
    end
    for (root, dirs, files) in walkdir(module_name)
        for file in files
            path = joinpath(root, file)
            if endswith(path, ".toml") && !endswith(dirname(path), ".tmp")
                check_results(path)
            end
        end
    end
end

using MarbleRun
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

        testname = first(splitext(basename(tomlfile)))
        if endswith(testname, "Error")
            E = Base.eval(MarbleRun, Symbol(last(split(testname, "_"))))
            @test_throws E MarbleRun.main(tomlfile)
            return
        else
            @time MarbleRun.main(tomlfile)
        end

        input = MarbleRun.readinputfile(tomlfile)
        for phase_index in 1:length(input.Phase)
            proj_dir = input.project
            output_dir = joinpath(input.Output.directory, string(phase_index))

            # for restart/readpointstate
            # remove suffix because the expected results are the same as the original simulation
            testname_ref = replace(replace(testname, "_restart"=>""), "_readpointstate"=>"") * "_$phase_index"

            # vtk points
            check_vtkpoints(
                output_dir,
                joinpath(proj_dir, "output", "$testname_ref.vtu"),
                0.05 * input.General.grid_space;
                fix_results
            )

            # snapshots file
            if input.Output.snapshots == true
                root, _, files = only(walkdir(joinpath(output_dir, "snapshots")))
                count = 0
                for file in sort(files, lt = natural)
                    check_snapshot(joinpath(root, "snapshot$count"); fix_results)
                    count += 1
                end
            end
            input.Output.snapshot_first == true && check_snapshot(joinpath(output_dir, "snapshots", "snapshot_first"); fix_results)
            input.Output.snapshot_last  == true && check_snapshot(joinpath(output_dir, "snapshots", "snapshot_last"); fix_results)

            # test history.csv if exists
            if isfile(joinpath(output_dir, "history.csv"))
                src = joinpath(output_dir, "history.csv")
                expected = joinpath(proj_dir, "output", "$testname_ref.csv")
                check_history(src, expected; fix_results)
            end

            # dirichlet
            if any(d->d.output, input.BoundaryCondition.Dirichlet)
                Dirichlet = input.BoundaryCondition.Dirichlet
                for i in eachindex(Dirichlet)
                    dirichlet = Dirichlet[i]
                    if dirichlet.output
                        src = joinpath(output_dir, "dirichlet", "$i", "history.csv")
                        expected = joinpath(proj_dir, "output", "$(testname_ref)_diriclet_$i.csv")
                        check_history(src, expected; fix_results)
                    end
                end
            end

            # rigidbodies
            if any(r->r.output, input.RigidBody) && startswith(tomlfile, "FreeRun")
                RigidBody = input.RigidBody
                for i in eachindex(RigidBody)
                    rigidbody = RigidBody[i]
                    if rigidbody.output
                        src = joinpath(output_dir, "rigidbodies", "$i", "history.csv")
                        expected = joinpath(proj_dir, "output", "$(testname_ref)_rigidbody_$i.csv")
                        check_history(src, expected; fix_results)
                    end
                end
            end
        end
    end
end

function check_vtkpoints(output_dir::String, expected::String, ϵ::Real; fix_results::Bool)
    # extract last vtk file
    vtk_file = joinpath(
        output_dir,
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
        cp(vtk_file, expected; force = true)
    else
        expected_points = get_points(VTKFile(expected))
        result_points = get_points(VTKFile(vtk_file))
        @test size(expected_points) == size(result_points)
        @test all(eachindex(expected_points)) do i
            norm(expected_points[i] - result_points[i]) < ϵ
        end
    end
end

function check_snapshot(file::String; fix_results::Bool)
    if fix_results
    else
        @test isfile(file)
        @test deserialize(file) isa NamedTuple
    end
end

function check_history(src::String, expected::String; fix_results::Bool)
    if fix_results
        cp(src, expected; force = true)
    else
        # check results
        history = CSV.File(src)
        output = CSV.File(expected) # expected output
        for name in propertynames(output)
            history_col = history[name]
            output_col = output[name][end-length(history_col)+1:end] # extract results for restart case
            @test output_col ≈ history_col atol=1e-4 rtol=0.01
        end
    end
end

@testset "Utilities" begin
    @testset "replace_version" begin
        # replace_version
        str = """
        # MarbleRun = "0.15"
        MarbleRun = "0.15" # should be replaced only this line
        """
        @test MarbleRun.replace_version(str, v"0.15.2") == """
        # MarbleRun = "0.15"
        MarbleRun = "0.15.2" # "0.15" is replaced by MarbleRun
        """
    end
end

@testset "$module_name" for module_name in ("GroundPenetration", "FreeRun",)
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

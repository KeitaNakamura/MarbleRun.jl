using MarbleRun
using Test

using LinearAlgebra # norm
using CSV           # history file
using Serialization # snapshots

# vtk
using ReadVTK
using NaturalSort

const FIX_RESULTS = false

const EXAMPLES = [
    "GroundPenetration" => [
        "OpenEndedPile/OpenEndedPile.toml",
        "OpenEndedPile/OpenEndedPile_restart.toml",
        "Spudcan/Spudcan.toml",
        "Spudcan/Spudcan_UndefKeyError.toml",
        "StripFooting/DruckerPrager.toml",
        "StripFooting/VonMises.toml",
    ],
    "FreeRun" => [
        "DamBreak/DamBreak.toml",
        "DamBreak/DamBreak_UnsupportedKeyError.toml",
        "SandColumn/SandColumn.toml",
        "SandColumn/SandColumn_readpointstate.toml",
        "SandColumn/SandColumn_restart.toml",
        "SandColumn/SandColumnWithSquareObject.toml",
        "SandColumn/SandColumnWithTriangleObject.toml",
        "SettlingDisk/SettlingDisk.toml",
        "StripFooting/StripFooting.toml",
    ],
]

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
                0.01 * input.General.grid_space
            )

            # snapshots file
            if input.Output.snapshots == true
                root, _, files = only(walkdir(joinpath(output_dir, "snapshots")))
                count = 0
                for file in sort(files, lt=natural)
                    check_snapshot(joinpath(root, "snapshot$count"))
                    count += 1
                end
            end
            input.Output.snapshot_first == true && check_snapshot(joinpath(output_dir, "snapshots", "snapshot_first"))
            input.Output.snapshot_last  == true && check_snapshot(joinpath(output_dir, "snapshots", "snapshot_last"))

            # test history.csv if exists
            if isfile(joinpath(output_dir, "history.csv"))
                results = joinpath(output_dir, "history.csv")
                expected = joinpath(proj_dir, "output", "$testname_ref.csv")
                check_history(results, expected)
            end

            # dirichlet
            if any(d->d.output, input.BoundaryCondition.Dirichlet)
                Dirichlet = input.BoundaryCondition.Dirichlet
                for i in eachindex(Dirichlet)
                    dirichlet = Dirichlet[i]
                    if dirichlet.output
                        results = joinpath(output_dir, "dirichlet", "$i", "history.csv")
                        expected = joinpath(proj_dir, "output", "$(testname_ref)_diriclet_$i.csv")
                        check_history(results, expected)
                    end
                end
            end

            # rigidbodies
            if any(r->r.output, input.RigidBody) && startswith(tomlfile, "FreeRun")
                RigidBody = input.RigidBody
                for i in eachindex(RigidBody)
                    rigidbody = RigidBody[i]
                    if rigidbody.output
                        results = joinpath(output_dir, "rigidbodies", "$i", "history.csv")
                        expected = joinpath(proj_dir, "output", "$(testname_ref)_rigidbody_$i.csv")
                        check_history(results, expected)
                    end
                end
            end
        end
    end
end

function check_vtkpoints(output_dir::String, expected::String, ϵ::Real)
    # extract last vtk file
    vtk_file = joinpath(
        output_dir,
        "paraview",
        sort(
            filter(
                file -> endswith(file, "_1.vtu"),
                only(walkdir(joinpath(output_dir, "paraview")))[3]
            ),
            lt=natural,
        )[end],
    )
    if FIX_RESULTS
        cp(vtk_file, expected; force=true)
    else
        expected_points = get_points(VTKFile(expected))
        result_points = get_points(VTKFile(vtk_file))
        @test size(expected_points) == size(result_points)
        @test all(eachindex(expected_points)) do i
            norm(expected_points[i] - result_points[i]) < ϵ
        end
    end
end

function check_snapshot(file::String)
    if FIX_RESULTS
    else
        @test isfile(file)
        @test deserialize(file) isa NamedTuple
    end
end

function check_history(results_csv::String, expected_csv::String)
    if FIX_RESULTS
        cp(results_csv, expected_csv; force=true)
    else
        # check results
        results = CSV.File(results_csv)
        expected = CSV.File(expected_csv)
        for name in propertynames(expected)
            results_col = results[name]
            expected_col = expected[name][end-length(results_col)+1:end] # extract results for restart case
            @test results_col ≈ expected_col atol=1e-8 rtol=0.01
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

# run examples
@testset "$mod" for (mod, list) in EXAMPLES
    for file in list
        path = joinpath(mod, file)
        rm(first(splitext(path)) * ".tmp"; force=true, recursive=true) # remove old folder
        check_results(path)
    end
end

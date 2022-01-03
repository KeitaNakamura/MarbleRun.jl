module Injection

using PoingrSimulator.GeometricObjects
using DelimitedFiles

const rigidbody_center_0 = Ref(Vec(NaN, NaN))

function main_output(args)
    history_file = joinpath(args.output_dir, "history.csv")

    grid = args.grid
    rigidbody = args.rigidbody

    if args.output_index == 0
        open(history_file, "w") do io
            writedlm(io, ["disp" "force"], ',')
        end
        rigidbody_center_0[] = centroid(rigidbody)
    end

    open(history_file, "a") do io
        disp = abs(centroid(rigidbody)[2] - rigidbody_center_0[][2])
        force = -sum(grid.state.fc)[2]
        writedlm(io, [disp force], ',')
    end
end

end

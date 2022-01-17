module Injection

using PoingrSimulator.GeometricObjects

const rigidbody_center_0 = Ref(Vec(NaN, NaN))

function main_output(args)
    history_file = joinpath(args.INPUT.Output.directory, "history.csv")

    grid = args.grid
    rigidbody = args.rigidbody

    if args.output_index == 0
        write(history_file, "disp,force\n")
        rigidbody_center_0[] = centroid(rigidbody)
    end

    open(history_file, "a") do io
        disp = abs(centroid(rigidbody)[2] - rigidbody_center_0[][2])
        force = -sum(grid.state.fc)[2]
        write(io, join([disp, force], ",") * "\n")
    end
end

end

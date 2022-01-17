module Injection

using PoingrSimulator.GeometricObjects

function main_output(args)
    history_file = joinpath(args.INPUT.Output.directory, "history.csv")

    t = args.t
    rigidbody = only(args.rigidbodies)

    if args.output_index == 0
        write(history_file, join(["t", "x", "y", "v_x", "v_y"], ",") * "\n")
    end

    open(history_file, "a") do io
        x, y = centroid(rigidbody)
        vx, vy = rigidbody.v
        open(history_file, "a") do io
            write(io, join([t, x, y, vx, vy], ",") * "\n")
        end
    end
end

end

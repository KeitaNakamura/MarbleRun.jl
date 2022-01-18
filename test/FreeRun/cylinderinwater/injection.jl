module Injection

using PoingrSimulator.GeometricObjects

function main_output_initialize(args)
    INPUT = args.INPUT

    history_file = joinpath(INPUT.Output.directory, "history.csv")
    write(history_file, join(["t", "x", "y", "v_x", "v_y"], ",") * "\n")
end

function main_output(args)
    INPUT = args.INPUT
    rigidbody = only(args.rigidbodies)
    t = args.t

    history_file = joinpath(INPUT.Output.directory, "history.csv")
    open(history_file, "a") do io
        x, y = centroid(rigidbody)
        vx, vy = rigidbody.v
        open(history_file, "a") do io
            write(io, join([t, x, y, vx, vy], ",") * "\n")
        end
    end
end

end

module Injection

using PoingrSimulator.GeometricObjects

function main_output_initialize(args)
    input = args.input

    history_file = joinpath(input.Output.directory, "history.csv")
    write(history_file, join(["t", "x", "y", "v_x", "v_y"], ",") * "\n")
end

function main_output(args)
    input = args.input
    rigidbody = only(args.rigidbodies)
    t = args.t

    history_file = joinpath(input.Output.directory, "history.csv")
    open(history_file, "a") do io
        x, y = centroid(rigidbody)
        vx, vy = rigidbody.v
        open(history_file, "a") do io
            write(io, join([t, x, y, vx, vy], ",") * "\n")
        end
    end
end

end

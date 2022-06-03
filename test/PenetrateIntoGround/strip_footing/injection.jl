module Injection

using MarbleBot.GeometricObjects

function main_output_initialize(args)
    input = args.input

    history_file = joinpath(input.Output.directory, "history.csv")
    write(history_file, "disp,force\n")
end

function main_output(args)
    input = args.input
    gridstate = args.gridstate
    rigidbody = args.rigidbody
    rigidbody0 = args.rigidbody0

    history_file = joinpath(input.Output.directory, "history.csv")
    open(history_file, "a") do io
        disp = abs(centroid(rigidbody)[2] - centroid(rigidbody0)[2])
        force = -sum(gridstate.fc)[2]
        write(io, join([disp, force], ",") * "\n")
    end
end

end

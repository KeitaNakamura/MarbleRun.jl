module PoingrSimulator

using Poingr
using GeometricObjects
using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

export PenetrateIntoGround

include("utils.jl")
include("PenetrateIntoGround.jl")

end # module

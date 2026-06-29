module MicroSwimmersPlots

using MicroSwimmers
using Makie
using LinearAlgebra
using GeometryBasics
using Parameters
using Meshing
using StaticArrays
using Interpolations

include("mesh_gen.jl")
include("viz.jl")
include("animations.jl")
include("exports.jl")

end # module
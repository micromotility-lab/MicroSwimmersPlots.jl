module MicroSwimmersPlots

using MicroSwimmers
using Makie
using LinearAlgebra
using GeometryBasics
using Parameters
using Meshing
using StaticArrays

include("mesh_gen.jl")
include("viz.jl")
include("animations.jl")
include("exports.jl")

set_theme!(theme_dark())

end # module
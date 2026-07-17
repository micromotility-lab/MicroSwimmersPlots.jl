module MicroSwimmersPlots

using MicroSwimmers
using Makie
using LinearAlgebra
using GeometryBasics
using Parameters
using Meshing
using StaticArrays
using Interpolations
using Accessors

include("implicit.jl")
include("mesh_gen.jl")
include("viz2.jl")
include("design2.jl")
include("animations.jl")
include("exports.jl")

end # module
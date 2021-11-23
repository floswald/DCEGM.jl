module DCEGM


using JSON
using Interpolations
using Distributions: Normal
using FastGaussQuadrature
using Roots
using Plots
using Colors
using Statistics
using Random
using Interact
using DataStructures: OrderedDict
using InvertedIndices: Not
using Blink
using BenchmarkTools


import Base: +, -, *, /, ==, isless, in
import Base: iterate, typemin, isapprox, convert
import Base.size, Base.getindex, Base.setindex!, Base.eltype
import Base.prepend!,
       Base.append!,
       Base.insert!,
       Base.delete!,
       Base.sort!,
       Base.length,
       Base.show,
       Base.IndexStyle


# Types
export MLine, Point, Envelope, yisless

# methods
export interp, splitat, getx, gety, gets, splitLine, getr

# includes for MultiLine
include("line.jl")
include("envelope.jl")

# # includes for dcegm
include("param.jl")
include("utils.jl")
include("minimal_EGM.jl")
include("dc_algo.jl")
include("simulate.jl")

include("plotting.jl")
include("bench.jl")
include("bk.jl")
include("interact.jl")


end # module

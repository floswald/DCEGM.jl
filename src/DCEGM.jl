module DCEGM


using MultiLine
using JSON
using Interpolations: interpolate, Gridded, Linear
using MiniLogging
using Distributions: Normal
using FastGaussQuadrature
using Plots

pyplot()

include("param.jl")
include("utils.jl")
include("dc_algo.jl")


end # module

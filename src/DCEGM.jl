module DCEGM


using MultiLine
using JSON
using Interpolations: interpolate, Gridded, Linear, extrapolate
using MiniLogging
using Distributions: Normal
using FastGaussQuadrature
using Plots

plotlyjs()

include("param.jl")
include("utils.jl")
include("dc_algo.jl")


end # module

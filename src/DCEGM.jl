module DCEGM


using JSON
using Interpolations: interpolate, Gridded, Linear, extrapolate
using Distributions: Normal
using FastGaussQuadrature
using Roots
using Plots
using ASTInterpreter2
gr()

# Types
export Line, Point, Envelope

# methods
export interp, splitat,upper_env!, getx, gety, gets, splitLine, getr

import Base.size, 
       Base.getindex, 
       Base.setindex!, 
       Base.eltype,
       Base.prepend!,
       Base.append!,
       Base.insert!,
       Base.delete!,
       Base.sort!,
       Base.length,
       Base.show

# includes for MultiLine
include("line.jl")
include("envelope.jl")
include("plotting.jl")

# includes for dcegm
include("param.jl")
include("utils.jl")
include("dc_algo.jl")


end # module

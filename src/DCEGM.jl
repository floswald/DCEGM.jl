module DCEGM


using JSON
using Interpolations: interpolate, Gridded, Linear, extrapolate
using Distributions: Normal
using FastGaussQuadrature
using Roots
using Plots

gr()

import Base: +, -, *, /, ==
import Base: start, next, done

# Types
export Line, Point, Envelope

# methods
export interp, splitat,upper_env!, getx, gety, gets, splitLine, getr

import Base.size, 
       Base.getindex, 
       Base.length, 
       Base.setindex!, 
       Base.eltype,
       Base.prepend!,
       Base.append!,
       Base.insert!,
       Base.delete!,
       Base.sort!,
       Base.show,
       Base.isless

# includes for MultiLine
include("line.jl")
# include("envelope.jl")

# # includes for dcegm
# include("param.jl")
# include("utils.jl")
# include("dc_algo.jl")
# include("plotting.jl")


function doi()

        x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
       L1 = Line(x,y)
       e = splitLine(L1)
       # upper_env!(e)
       # plot(e)
       e
end
end # module

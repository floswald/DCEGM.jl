module DCEGM


using JSON
using Interpolations
using Distributions: Normal
using FastGaussQuadrature
using Roots
using Plots
using Colors
using Statistics
gr()


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
include("dc_algo.jl")
include("simulate.jl")

include("plotting.jl")
include("bench.jl")


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

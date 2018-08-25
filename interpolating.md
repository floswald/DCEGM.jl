# Interpolating `Point`s

I made this doc with the command
```julia
using Literate; Literate.markdown("examples/interpolating.jl",".";documenter=false)
```
This document shows how to interpolate a `Vector{Point{T}}` instead of two standard `Vector{T}`:

```julia
module Points
using Interpolations
using BenchmarkTools

import Base: +, -, *, /

struct Point{T}
	x::T
	y::T
end
```

Having defined the `Point` struct, we next define some arithmetic operations that will allow us to use the type with the `Interpolations.jl` package:

```julia
(+)(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y)
(-)(p1::Point, p2::Point) = Point(p1.x-p2.x, p1.y-p2.y)
(*)(n::Number, p::Point) = Point(n*p.x, n*p.y)
(*)(p::Point, n::Number) = n*p
(/)(p::Point, n::Number) = Point(p.x/n, p.y/n)

Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.promote_rule(::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{T1}, ::Type{Point{T2}}) where {T1<:Number,T2} = Point{promote_type(T1,T2)}

# (+)(p1::Point, p2::Point) = p1.y+p2.y
# (-)(p1::Point, p2::Point) = p1.y-p2.y
# (*)(n::Number, p::Point) = n*p.y
# (*)(p::Point, n::Number) = n*p
# (/)(p::Point, n::Number) =  p.y/n

# @noinline checksame(p1, p2) = p1.x == p2.x || throw(ArgumentError("arithmetic supported only for constant x, got $p1 and $p2"))
# function (+)(p1::Point, p2::Point)
#     # checksame(p1, p2)
#     Point(p1.x, p1.y + p2.y)
# end
# function (-)(p1::Point, p2::Point)
#     # checksame(p1, p2)
#     Point(p1.x, p1.y - p2.y)
# end
# (*)(n::Number, p::Point) = Point(n*p.x, n*p.y)
# (/)(n::Number, p::Point) = Point(n/p.x, n/p.y)
# (*)(p::Point, n::Number) = n*p
```

Then we can move on to test the performance of this new way of interpolating. Notice that there are two concerns in the speed comparison:
1. how long does it take to interpolate a `Vector{Point{T}}` vs 2 `Vector{T}`s?
1. how long does it take to convert data in `Vector{Point{T}}` into 2 `Vector{T}`s?

```julia
n = 1_000_000
m = 10_000
```

cast data as vector of Point

```julia
data = reinterpret(Point{Float64},hcat(collect(1:n),rand(n))',(n,))
# data to test interpolation on
rdata = rand(m)
```

convert vec of points to separate x and y vectors

```julia
x = [i.x for i in data]
y = [i.y for i in data]
```

my favourite version
but also interpolates the x-grid :-(

```julia
function itp_points(d,rd)
	# itp = scale(interpolate(d, BSpline(Quadratic(Flat())),OnGrid()),linspace(d[1].x,d[end].x,length(d)))
	itp = scale(interpolate(d, BSpline(Linear()),OnGrid()),linspace(d[1].x,d[end].x,length(d)))
	out = similar(d)
	for i in 1:length(rd)
		out[i] = itp[rd[i]]
	end
	out
end
```

traditional version

```julia
# only interpolate the y vector
function itp_x_y(x,y,rd)
	itp = scale(interpolate(y, BSpline(Linear()),OnGrid()),linspace(x[1],x[end],length(x)))
	out = similar(rd)
	for i in 1:length(rd)
		out[i] = itp[rd[i]]
	end
	out
end
```

actual cost of traditional:
need to first convert Point to vector!
an in fact that conversion would have to happen FOR EACH interpolation
in the below loop... so much more costly in the end!

```julia
function itp_pointsxy(d,rd)
	x = [i.x for i in data]
	y = [i.y for i in data]
	itp = scale(interpolate(y, BSpline(Linear()),OnGrid()),linspace(x[1],x[end],length(x)))
	out = similar(rd)
	for i in 1:length(rd)
		out[i] = itp[rd[i]]
	end
	out
end
```

compile

```julia
itp_points(data[1:10],rdata[1:5])
itp_x_y(x[1:10],y[1:10],rdata[1:5])
itp_pointsxy(data[1:10],rdata[1:5])
```

BenchmarkTools

```julia
t1 = @benchmark itp_points($data,$rdata)
t2 = @benchmark itp_x_y($x,$y,$rdata)
t3 = @benchmark itp_pointsxy($data,$rdata)

1+2

end

# We see here the timings for t1 and t2:
```

```julia
julia v0.6.2> Points.t1
BenchmarkTools.Trial:
  memory estimate:  30.52 MiB
  allocs estimate:  8
  --------------
  minimum time:     6.088 ms (9.38% GC)
  median time:      8.402 ms (34.32% GC)
  mean time:        8.675 ms (34.34% GC)
  maximum time:     26.954 ms (28.22% GC)
  --------------
  samples:          576
  evals/sample:     1

julia v0.6.2> Points.t2
BenchmarkTools.Trial:
  memory estimate:  7.71 MiB
  allocs estimate:  8
  --------------
  minimum time:     2.007 ms (0.00% GC)
  median time:      2.213 ms (0.00% GC)
  mean time:        3.514 ms (37.21% GC)
  maximum time:     13.390 ms (77.77% GC)
  --------------
  samples:          1419
  evals/sample:     1
  ```

# a hacked arithmetic for this type

```julia
# now lets redefine the module, in particular let's change the arithmetics:

module Points
using Interpolations
using BenchmarkTools

import Base: +, -, *, /

struct Point{T}
	x::T
	y::T
end

# (+)(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y)
# (-)(p1::Point, p2::Point) = Point(p1.x-p2.x, p1.y-p2.y)
# (*)(n::Number, p::Point) = Point(n*p.x, n*p.y)
# (*)(p::Point, n::Number) = n*p
# (/)(p::Point, n::Number) = Point(p.x/n, p.y/n)

Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.promote_rule(::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{T1}, ::Type{Point{T2}}) where {T1<:Number,T2} = Point{promote_type(T1,T2)}

(+)(p1::Point, p2::Point) = p1.y+p2.y
(-)(p1::Point, p2::Point) = p1.y-p2.y
(*)(n::Number, p::Point) = n*p.y
(*)(p::Point, n::Number) = n*p
(/)(p::Point, n::Number) =  p.y/n

n = 1_000_000
m = 10_000
data = reinterpret(Point{Float64},hcat(collect(1:n),rand(n))',(n,))
rdata = rand(m)
x = [i.x for i in data]
y = [i.y for i in data]

function itp_points(d,rd)
	# itp = scale(interpolate(d, BSpline(Quadratic(Flat())),OnGrid()),linspace(d[1].x,d[end].x,length(d)))
	itp = scale(interpolate(d, BSpline(Linear()),OnGrid()),linspace(d[1].x,d[end].x,length(d)))
	out = similar(rd)
	for i in 1:length(rd)
		out[i] = itp[rd[i]]
	end
	out
end

function itp_x_y(x,y,rd)
	itp = scale(interpolate(y, BSpline(Linear()),OnGrid()),linspace(x[1],x[end],length(x)))
	out = similar(rd)
	for i in 1:length(rd)
		out[i] = itp[rd[i]]
	end
	out
end

function itp_pointsxy(d,rd)
	x = [i.x for i in data]
	y = [i.y for i in data]
	itp = scale(interpolate(y, BSpline(Linear()),OnGrid()),linspace(x[1],x[end],length(x)))
	out = similar(rd)
	for i in 1:length(rd)
		out[i] = itp[rd[i]]
	end
	out
end

itp_points(data[1:10],rdata[1:5])
itp_x_y(x[1:10],y[1:10],rdata[1:5])
itp_pointsxy(data[1:10],rdata[1:5])

t1 = @benchmark itp_points($data,$rdata)
t2 = @benchmark itp_x_y($x,$y,$rdata)
t3 = @benchmark itp_pointsxy($data,$rdata)

end
```

This shows that the *hacked* method is a bit faster, and has much less memory allocs.

```julia
julia v0.6.2> Points.t1
BenchmarkTools.Trial:
  memory estimate:  15.34 MiB
  allocs estimate:  8
  --------------
  minimum time:     3.760 ms (0.00% GC)
  median time:      6.865 ms (13.17% GC)
  mean time:        5.885 ms (26.24% GC)
  maximum time:     20.627 ms (32.43% GC)
  --------------
  samples:          848
  evals/sample:     1

julia v0.6.2> Points.itp_points(Points.data[1:10],Points.rdata[1:5])
5-element Array{Float64,1}:
 0.178303
 0.171897
 0.184147
 0.179542
 0.159741

julia v0.6.2> Points.itp_x_y(Points.x[1:10],Points.y[1:10],Points.rdata[1:5])
5-element Array{Float64,1}:
 0.178303
 0.171897
 0.184147
 0.179542
 0.159741
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


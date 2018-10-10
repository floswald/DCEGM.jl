
struct Point{T}
    x::T
    y::T
end
eltype(p::Point) = eltype(p.x) 

function show(io::IO, ::MIME"text/plain", p::Point{T}) where T
    print(io,"Point of type $T:\n")
    print(io,"      x = $(p.x)\n")
    print(io,"      y = $(p.y)\n")
end

show(io::IO,p::Point{T}) where T = print(io,"($(p.x),$(p.y)")

# arithmetics
(+)(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y)
(-)(p1::Point, p2::Point) = Point(p1.x-p2.x, p1.y-p2.y)
(*)(n::Number, p::Point) = Point(n*p.x, n*p.y)
(*)(p::Point, n::Number) = n*p
(/)(p::Point, n::Number) = Point(p.x/n, p.y/n)
(==)(p1::Point, p2::Point) = (p1.x == p2.x) && (p1.y == p2.y)

"""
# Line

A `Line` is a vector of `Point`s. The x-coordinates of the points can be irregularly spaced.

## Fields

* `v`: Vector of `Point`s
* `n`: number of points in line
* `xvec`: a vector of the x values for gridded interpolation
* `xrange`: the xrange of the line, i.e. the range of `x`.
* `yrange`: the yrange of the line, i.e. the range of `x`.
"""
mutable struct Line{T<:Number} <: AbstractArray{T<:Number,1}
    v::Vector{Point{T}}
    n::Int
    xvec ::Vector{T}
    xrange::Tuple
    yrange::Tuple
    function Line(v::Vector{Point{T}}) where {T<:Number}
        this = new{T}()
        this.v = v
        this.n = length(v)
        this.xvec = [i.x for i in v]
        y = [i.y for i in v]
        this.xrange = extrema(this.xvec)
        this.yrange = extrema(y)
        return this
    end
    function Line(x::Vector{T},y::Vector{T}) where {T<:Number}
        this = new{T}()
        n = length(x)
        @assert n == length(y)
        this.xvec = copy(x)
        this.v = [Point(x[i],y[i]) for i in 1:n]
        this.n = length(this.v)
        this.xrange = extrema(this.xvec)
        this.yrange = extrema(y) 
        return this
    end
end

# printing
function show(io::IO, ::MIME"text/plain", L::Line{T}) where {T<:Number}
    print(io,"$T Line\n")
    print(io,"number of points: $(L.n)\n")
    print(io,"xrange: $(L.xrange)\n")
    print(io,"yrange: $(L.yrange)\n")
end
show(io::IO,L::Line{T}) where {T<:Number} = print(io,"$(L.n) point $T Line")

# indexing
eltype(l::Line) = eltype(l.v) 
size(l::Line) = (l.n,)
length(l::Line) = l.n
function getindex(l::Line,i)
    @boundscheck checkbounds(l,i)
    l.v[i...]
end
function setindex!(l::Line{T},v::Point{T},i::Int) where {T<:Number}
    l[i] = v
    m.xvec = [i.x for i in m.v]
end
endof(l::Line) = l.n


# iteration
Base.start(::Line) = 1
Base.next(L::Line,state) = (L[state],state+1)
Base.done(L::Line,state) = state > length(L)

# min max
function isless(p1::Point{T},p2::Point{T})  where T
    p1.x < p2.x
end
function max_x(l::Line{T}) where {T<:Number}
    lo = typemin(T)
    for i in l
        if i.x > lo
            lo = i.x
        end
    end
    return lo
end
function min_x(l::Line{T}) where {T<:Number}
    hi = typemax(T)
    for i in l
        if i.x < hi
            hi = i.x
        end
    end
    return hi
end
function max_y(l::Line{T}) where {T<:Number}
    lo = typemin(T)
    for i in l
        if i.y > lo
            lo = i.y
        end
    end
    return lo
end
function min_y(l::Line{T}) where {T<:Number}
    hi = typemax(T)
    for i in l
        if i.y < hi
            hi = i.y
        end
    end
    return hi
end

# function getindex(l::Line,i::UnitRange{Int})
#     Line(l.x[i],l.y[i])
# end
# function setindex!(l::Line{T},x::T,y::T,i::Int) where {T<:Number}
#     l.x[i] = x
#     l.y[i] = y
# end
# function setindex!(l::Line{T},x::Vector{T},y::Vector{T},i::UnitRange{Int}) where {T<:Number}
#     l.x[i] = x
#     l.y[i] = y
# end

# interpolating a line

"""
   interp(l::Line{T},ix::Vector{T},extrap::Bool=false) where {T<:Number}

Interpolate a `Line` on a vector of valuels `x`
"""
function interp(l::Line{T},ix::Vector{T},extrap::Bool=false) where {T<:Number}
    # whenever 
    xex = extrema(ix)
    # @debug(logger,"interpolating $ix ")
    if xex[1] < l.xrange[1] || xex[2] > l.xrange[2]
        if extrap 
            # itp = scale(interpolate(l.v, BSpline(Linear()),OnGrid()),linspace(d[1].x,d[end].x,length(d)))
            itp = extrapolate(interpolate((l.xvec,),l.v,Gridded(Linear())),Linear())
        else
            itp = extrapolate(interpolate((l.xvec,),l.v,Gridded(Linear())),-Inf)
        end
    else
        itp = interpolate((l.xvec,),l.v,Gridded(Linear()))
    end
    return itp[ix]
end 
function interp(e::Array{Line{T}},ix::Vector{T};extrap::Bool=false) where {T<:Number}
    y = zeros(T,length(e),length(ix))
    for i in eachindex(e)
        y[i,:] = interp(e[i],ix,extrap)
    end
    return y
end 

# appending, prepending , deleting and splitting at

"prepend `Point`s to a Line"
function prepend!(m::Line{T},p::Vector{Point{T}}) where T
    prepend!(m.v,p)
    reconfigure!(m)
end

function reconfigure!(m::Line)
    # after having updated some objects, need to recompute n
    m.xvec = [i.x for i in m.v]
    m.n = length(m.v)
    m.xrange = extrema(m.xvec)
    m.yrange = (min_y(m),max_y(m))
end

"delete an index"
function delete!(m::Line,idx)
    deleteat!(m.v,idx)
    reconfigure!(m)
end

"append points to a Line"
function append!(m::Line{T},p::Vector{Point{T}}) where T
    append!(m.v,p)
    reconfigure!(m)
end

"insert a single value at an interior index"
function insert!(m::Line{T},v::Point{T},idx::Int) where T
    insert!(m.v,idx,v)
    reconfigure!(m)
end

"""
    splitat(m::Line,j::Int,repeat_boundary::Bool=true)

Splits a `Line` object after given index and returns 2 new `Line`s as a tuple. If `repeat_boundary` is true, then the separating index is the first point of the second new `Line`.
"""
function splitat(m::Line,j::Int,repeat_boundary::Bool=true)
    m1 = Line(m.v[1:j])
    if repeat_boundary
        m2 = Line(m.v[j:end])
    else
        m2 = Line(m.v[j+1:end])
    end
    return (m1,m2)
end

"sort a `Line` along x-grid"
function sortx!(m::Line)
    sort!(m.v)  # sorts v
    reconfigure!(m)
end


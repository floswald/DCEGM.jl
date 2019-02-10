
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

show(io::IO,p::Point{T}) where T = print(io,"($(p.x),$(p.y))")

# arithmetics
(+)(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y)
(-)(p1::Point, p2::Point) = Point(p1.x-p2.x, p1.y-p2.y)
(*)(n::Number, p::Point) = Point(n*p.x, n*p.y)
(*)(p::Point, n::Number) = n*p
(/)(p::Point, n::Number) = Point(p.x/n, p.y/n)


# comparison
(==)(p1::Point, p2::Point) = (p1.x == p2.x) && (p1.y == p2.y)
isapprox(p1::Point, p2::Point) = isapprox(p1.x, p2.x) && isapprox(p1.y, p2.y)
isless(p1::Point,p2::Point) = p1.x < p2.x   # caution we sort on x only here! 


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

# typemin
eltype(l::Line) = eltype(l.v) 
function typemin(L::Line{T}) where {T<:Number}
    typemin(eltype(L))
end

# iteration
Base.start(::Line) = 1
Base.next(L::Line,state) = (L[state],state+1)
Base.done(L::Line,state) = state > length(L)

# printing
function show(io::IO, ::MIME"text/plain", L::Line{T}) where {T<:Number}
    print(io,"$T Line\n")
    print(io,"number of points: $(L.n)\n")
    print(io,"xrange: $(L.xrange)\n")
    print(io,"yrange: $(L.yrange)\n")
end
show(io::IO,L::Line{T}) where {T<:Number} = print(io,"$(L.n) point $T Line")

getx(l::Line{T}) where {T<:Number}= T[l.v[i].x for i in 1:l.n] 
gety(l::Line{T}) where {T<:Number} = T[l.v[i].y for i in 1:l.n]

# indexing
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
function findmax_y(l::Line{T}) where {T<:Number}
    lo = typemin(T)
    ix = 0
    imax = 0
    for i in l
        ix += 1
        if i.y > lo
            lo = i.y
            imax = ix
        end
    end
    return (lo,imax)
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

Interpolate a `Line` on a vector of values `x`.
Importantly, this returns a new vector of `Point` (i.e. tuples of (x,y)).
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
            itp = extrapolate(interpolate((l.xvec,),l.v,Gridded(Linear())),Flat())
        end
    else
        itp = interpolate((l.xvec,),l.v,Gridded(Linear()))
    end
    return itp(ix)
end 

function interp(e::Array{Line{T}},ix::Vector{T};extrap::Bool=false) where {T<:Number}
    [Line(interp(e[i],ix,extrap)) for i in eachindex(e)]
end 

function linemax(e::Array{Line{T}}) where {T<:Number}
    out = zeros(Int,length(e[1]))
    rows = length(e)
    for j in eachindex(e[1])
        ic = 0
        iv = 0
        v = typemin(T)
        for row in 1:rows
            ic += 1
            if e[row][j].y > v
                v = e[row][j].y
                iv = ic
            end
        end
        out[j] = iv
    end
    return out
end 

function t1()
    n = 15
    x1 = collect(linspace(0,10,n))
    x2 = collect(linspace(-1,9,n))
    L1 = DCEGM.Line(x1,x1)
    L2 = Line(x2,ones(n)*5)
    e = DCEGM.Envelope([L1,L2])
    upper_env!(e)
    # ix = [1.1,1.2,6.0]
    # y=DCEGM.interp(e.L,ix)
    # linemax(y)
    return e
end

function t2()

    n = 15
    x1 = collect(linspace(0,10,n))
    x2 = collect(linspace(-1,9,n))
    x3 = collect(linspace(6.1,10,n))
    L1 = Line(x1,x1)
    L2 = Line(x2,ones(n)*5)
    L3 = Line(x3,(x3)*1.9 .- 7.0)
    e = Envelope([L1,L2,L3])
    upper_env!(e)
    removed!(e)
    e
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


"""
    intersect(L1::Line,L2::Line,s::Int)

Intersecting two lines returns a [`Point`](@ref)
"""
function intersect(L1::Line,L2::Line,s::Int)
    xlo,vlo = (L1[s].x,L1[s].y)
    xhi,vhi = (L2[s+1].x,L2[s+1].y)
    f_closure(z) = interp(L1,[z])[1].y - interp(L2,[z])[1].y
    if f_closure(xlo) * f_closure(xhi) > 0
        # not opposite signs, no zero to be found
        x_x = xlo
        v_x = vhi
        return Point(x_x,v_x)
    else
        x_x = fzero(f_closure,xlo,xhi)
        v_x = interp(L1,[x_x])[1].y
        return Point(x_x,v_x)
    end
end

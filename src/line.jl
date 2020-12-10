struct Point{T}
    x::T
    y::T
end
function convert(::Type{Point},x::Vector,y::Vector) where T
    @assert length(x) == length(y)
    [Point(x[ix],y[ix]) for ix in 1:length(x)]
end
eltype(p::Point) = eltype(p.x)

"get x coordinates from a vector of `Point`"
getx(v::Vector{Point{T}})  where T = T[v[i].x for i in 1:length(v)]

"get y coordinates from a vector of `Point`"
gety(v::Vector{Point{T}})  where T = T[v[i].y for i in 1:length(v)]

"get x-y coordinates from a vector of `Point`"
coords(v::Vector{Point{T}}) where T = [getx(v) gety(v)]

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
(==)(x::T, p2::Point{T}) where T = (x == p2.x)  # is x-coord part of this point?
isapprox(p1::Point, p2::Point) = isapprox(p1.x, p2.x) && isapprox(p1.y, p2.y)
isless(p1::Point,p2::Point) = p1.x < p2.x   # caution we sort on x only here!
in(x::T, p2::Vector{Point{T}}) where T = any(x .== p2)


# Promotion
Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.promote_rule(::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{T1}, ::Type{Point{T2}}) where {T1<:Number,T2} = Point{promote_type(T1,T2)}



"""
# MLine

A `MLine` is a vector of `Point`[@ref]. The x-coordinates of the points can be irregularly spaced.

## Fields

* `v`: Vector of `Point`s
* `iextrap`: indices of points that have been extrapolated
"""
mutable struct MLine{T} <: AbstractArray{T,1}
    v::Vector{Point{T}}
    iextrap::Vector{Int}
    function MLine(v::Vector{Point{T}}; extrapolated=nothing) where {T<:Number}
        this = new{T}()
        this.v = v
        if isnothing(extrapolated)
            this.iextrap = Int[]
        else
            this.iextrap = extrapolated
        end
        return this
    end
    function MLine(x::Vector{T},y::Vector{T}; extrapolated=nothing) where {T<:Number}
        this = new{T}()
        n = length(x)
        @assert n == length(y)
        this.v = [Point(x[i],y[i]) for i in 1:n]
        if isnothing(extrapolated)
            this.iextrap = Int[]
        else
            this.iextrap = extrapolated
        end
        return this
    end
end

# typemin
Base.eltype(l::MLine) = eltype(l.v)
function typemin(L::MLine{T}) where {T<:Number}
    typemin(eltype(L))
end

# indexing
Base.size(l::MLine) = (length(l.v),)
Base.length(l::MLine) = length(l.v)

function getindex(l::MLine,i::Int)
    # @boundscheck checkbounds(l,i)
    l.v[i]
end
Base.IndexStyle(::Type{<:MLine}) = IndexLinear()
getindex(l::MLine, I...) = l.v[I]
endof(l::MLine) = length(l.v)

# iteration
Base.iterate(L::MLine, state = 1) = state > length(L.v) ? nothing : (L[state], state + 1)

# printing
function show(io::IO, ::MIME"text/plain", L::MLine{T}) where {T<:Number}
    xvec = [i.x for i in L.v]
    print(io,"$T MLine\n")
    print(io,"number of points: $(length(L.v))\n")
    print(io,"xrange: $(round.(extrema(getx(L)),digits = 2))\n")
    print(io,"yrange: $(round.(extrema(gety(L)),digits = 2))\n")
end
show(io::IO,L::MLine{T}) where {T<:Number} = print(io,"$(length(L.v)) point $T MLine")

getv(l::MLine{T}) where T = l.v
getx(l::MLine{T}) where T = getx(l.v)
gety(l::MLine{T}) where T = gety(l.v)
getex(l::MLine{T}) where T = l.iextrap

coords(l::MLine{T}) where T = coords(l.v)
function unique!(l::MLine{T}) where T
    l.v = unique(l.v)
end

function floory!(l::MLine{T},yy::T) where {T<:Number}
    ix = findall(gety(l) .< yy)
    if length(ix) > 0
        newp = l.v[ix]
        for i in 1:length(newp)
            splice!(l.v,ix[i],[Point(newp[i].x,yy)])
        end
    end
end

# Array of MLine
function floory!(L::Array{MLine{T}},x::T) where {T<:Number}
    for i in eachindex(L)
        floory!(L[i],x)
    end
end
function gety(L::Vector{MLine{T}}) where {T<:Number}
    out = fill(zero(T),(length(L),length(L[1].v)))
    for i in eachindex(L)
        out[i,:] = gety(L[i])
    end
    out
end







function max_x(l::MLine{T}) where {T<:Number}
    lo = typemin(T)
    for i in l
        if i.x > lo
            lo = i.x
        end
    end
    return lo
end
function min_x(l::MLine{T}) where {T<:Number}
    hi = typemax(T)
    for i in l
        if i.x < hi
            hi = i.x
        end
    end
    return hi
end
function max_y(l::MLine{T}) where {T<:Number}
    lo = typemin(T)
    for i in l
        if i.y > lo
            lo = i.y
        end
    end
    return lo
end
function findmax_y(l::MLine{T}) where {T<:Number}
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
function min_y(l::MLine{T}) where {T<:Number}
    hi = typemax(T)
    for i in l
        if i.y < hi
            hi = i.y
        end
    end
    return hi
end


# interpolating a MLine

"""
   interp(l::MLine{T},ix::Vector{T};Bool::extrap = false) where {T<:Number}

Interpolate a `MLine` on a vector of values `x`.
Importantly, this returns a new vector of `Point` (i.e. tuples of (x,y), and not just a function value y).
"""
function interp(l::MLine{T},ix::Vector{T};extrap::Bool = true) where {T<:Number}
    # # whenever
    # xex = extrema(ix)
    # @debug(logger,"interpolating $ix ")
    xvec = getx(l)
    if !issorted(xvec)
        # println(xvec)
    end
    xrange = extrema(xvec)

    fi = findall((ix .< xrange[1]) .| (ix .> xrange[2]))

    # by default extrapolate all lines
    # but record which points have been extrapolated
    if length(fi) > 0
        if extrap
            itp = extrapolate(interpolate((xvec,),l.v,Gridded(Linear())),Line())
            out = MLine(itp(ix),extrapolated = fi)
        else
            out = MLine(zeros(T,length(ix)),zeros(T,length(ix)))
            for xi in fi
                out.v[xi] = Point(ix[xi],typemin(T))
            end
            for xi in setdiff(1:length(ix),fi)
            # for xi in collect((1:length(xvec)))[Not(fi)]
                out.v[xi] = interpolate((xvec,),l.v,Gridded(Linear()))(ix[xi])
            end
        end
    else
        itp = interpolate((xvec,),l.v,Gridded(Linear()))
        out = MLine(itp(ix))
    end

    return out
end

function interp(e::Array{MLine{T}},ix::Vector{T};extrap::Bool = true) where {T<:Number}
    [interp(e[i],ix,extrap=extrap) for i in eachindex(e)]
end


"""
    linemax(e::Array{MLine{T}}) where {T<:Number}

For an array of `MLine`s on identical support `xx`,
computes the index in `e` of the `MLine` where the `y`-value is highest for each `xx`.
One can imagine `e` as a matrix where each row represents the `y`-values from a
different `MLine`; this function returns the index of the column-wise maximum.

"""
function linemax(e::Array{MLine{T}}) where {T<:Number}
    out = zeros(Int,length(e[1]))
    rows = length(e)  # how many lines
    for j in eachindex(e[1])   # across columns
        ic = 0  # current index
        iv = 0  # index of best value
        v = typemin(T)
        for row in 1:rows
            ic += 1
            # if (e[row][j].y > v) && ( !(j ∈ e[row].iextrap) )   # if best value and not extrapolated
            if (e[row][j].y > v)    # if best value
                v = e[row][j].y
                iv = ic
            end
        end
        out[j] = iv
    end
    return out
end
# appending, prepending , deleting and splitting at

"prepend `Point`s to a `MLine`[@ref]"
function prepend!(m::MLine{T},p::Vector{Point{T}}) where T
    prepend!(m.v,p)
    # reconfigure!(m)
end


"delete an index"
function delete!(m::MLine,idx)
    deleteat!(m.v,idx)
    # reconfigure!(m)
end

"append points to a MLine"
function append!(m::MLine{T},p::Vector{Point{T}}) where T
    append!(m.v,p)
    # reconfigure!(m)
end

"insert a single value at an interior index"
function insert!(m::MLine{T},v::Point{T},idx::Int) where T
    insert!(m.v,idx,v)
    # reconfigure!(m)
end

"""
    splitat(m::MLine,j::Int,repeat_boundary::Bool=true)

Splits a `MLine` object after given index and returns 2 new `MLine`s as a tuple.
If `repeat_boundary` is true, then the separating index is the first point of the second new `MLine`.
Propagates the `iextrap` vector of indices to the newly formed `MLine`s.
"""
function splitat(m::MLine,j::Int,repeat_boundary::Bool=true)
    iex = getex(m)
    if length(iex) > 0
        # m had any extrapolated points
        # keep track of those in the new object
        m1 = MLine(m.v[1:j],extrapolated = any(iex .< j+1) ? iex[iex .< j+1] : nothing)

        if repeat_boundary
            m2 = MLine(m.v[j:end],extrapolated = any(iex .>= j) ? iex[iex .>= j] .- (j-1) : nothing)

        else
            m2 = MLine(m.v[j+1:end],extrapolated = any(iex .> j) ? iex[iex .> j] .- j : nothing)
        end
    else
        m1 = MLine(m.v[1:j])
        if repeat_boundary
            m2 = MLine(m.v[j:end])
        else
            m2 = MLine(m.v[j+1:end])
        end
    end
    return (m1,m2)
end

"sort a `MLine` along x-grid"
function sortx!(m::MLine)
    sort!(m.v)  # sorts v
    # reconfigure!(m)
end

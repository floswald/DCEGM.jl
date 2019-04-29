struct Point{T}
    x::T
    y::T
end
eltype(p::Point) = eltype(p.x) 
zero

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

# Promotion
Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.promote_rule(::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{Point{T1}}, ::Type{T2}) where {T1,T2<:Number} = Point{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{T1}, ::Type{Point{T2}}) where {T1<:Number,T2} = Point{promote_type(T1,T2)}



"""
# MLine

A `MLine` is a vector of `Point`s. The x-coordinates of the points can be irregularly spaced.

## Fields

* `v`: Vector of `Point`s
* `n`: number of points in line
* `xvec`: a vector of the x values for gridded interpolation
* `xrange`: the xrange of the line, i.e. the range of `x`.
* `yrange`: the yrange of the line, i.e. the range of `x`.
* `extrap`: whether this `MLine` should be extrapolated beyond it's domain. `true` by default
"""
mutable struct MLine{T} <: AbstractArray{T,1}
    v::Vector{Point{T}}
    n::Int
    xvec ::Vector{T}
    xrange::Tuple
    yrange::Tuple
    extrap::Bool
    function MLine(v::Vector{Point{T}}; extrap::Bool=true) where {T<:Number}
        this = new{T}()
        this.v = v
        this.n = length(v)
        this.xvec = [i.x for i in v]
        y = [i.y for i in v]
        this.xrange = extrema(this.xvec)
        this.yrange = extrema(y)
        this.extrap = extrap
        return this
    end
    function MLine(x::Vector{T},y::Vector{T}; extrap::Bool=true) where {T<:Number}
        this = new{T}()
        n = length(x)
        @assert n == length(y)
        this.xvec = copy(x)
        this.v = [Point(x[i],y[i]) for i in 1:n]
        this.n = length(this.v)
        this.xrange = extrema(this.xvec)
        this.yrange = extrema(y) 
        this.extrap = extrap
        return this
    end
end

# typemin
Base.eltype(l::MLine) = eltype(l.v) 
function typemin(L::MLine{T}) where {T<:Number}
    typemin(eltype(L))
end

# indexing
Base.size(l::MLine) = (l.n,)
Base.length(l::MLine) = l.n

function getindex(l::MLine,i)
    @boundscheck checkbounds(l,i)
    l.v[i...]
end
function setindex!(l::MLine{T},v::Point{T},i::Int) where {T<:Number}
    l[i] = v
    # reconfigure xvec
    m.xvec = getx(l)
end
endof(l::MLine) = l.n

# iteration
Base.iterate(L::MLine, state = 1) = state > L.n ? nothing : (L[state], state + 1)

# printing
function show(io::IO, ::MIME"text/plain", L::MLine{T}) where {T<:Number}
    print(io,"$T MLine\n")
    print(io,"number of points: $(L.n)\n")
    print(io,"xrange: $(L.xrange)\n")
    print(io,"yrange: $(L.yrange)\n")
end
show(io::IO,L::MLine{T}) where {T<:Number} = print(io,"$(L.n) point $T MLine")

getx(l::MLine{T}) where {T<:Number}= T[l.v[i].x for i in 1:l.n] 
gety(l::MLine{T}) where {T<:Number} = T[l.v[i].y for i in 1:l.n]




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
   interp(l::MLine{T},ix::Vector{T}) where {T<:Number}

Interpolate a `MLine` on a vector of values `x`.
Importantly, this returns a new vector of `Point` (i.e. tuples of (x,y)).
"""
function interp(l::MLine{T},ix::Vector{T}) where {T<:Number}
    # whenever 
    xex = extrema(ix)
    # @debug(logger,"interpolating $ix ")
    if l.extrap
        itp = extrapolate(interpolate((l.xvec,),l.v,Gridded(Linear())),Linear())
        out = MLine(itp(ix))
    else
        fi = findall((ix .< l.xrange[1]) .| (ix .> l.xrange[2]))

        # manually set out of bounds x to -Inf
        if length(fi) > 0
            out = MLine(zeros(T,length(ix)),zeros(T,length(ix)))
            for xi in fi
                out.v[xi] = Point(ix[xi],typemin(T))
            end
            for xi in setdiff(1:length(ix),fi)
                out.v[xi] = interpolate((l.xvec,),l.v,Gridded(Linear()))(ix[xi])
            end

        else
            itp = interpolate((l.xvec,),l.v,Gridded(Linear()))
            out = MLine(itp(ix))
        end
    end

    return out
end 

function interp(e::Array{MLine{T}},ix::Vector{T}) where {T<:Number}
    [interp(e[i],ix) for i in eachindex(e)]
end 
# function interp!(o::Matrix{Point{T}},e::Array{MLine{T}},ix::Vector{T};extrap::Bool=false) where {T<:Number}
#     for i in eachindex(e)
#         o[i,:] = MLine(interp(e[i],ix,extrap))
#     end
#     # [MLine(interp(e[i],ix,extrap)) for i in eachindex(e)]
# end 

function linemax(e::Array{MLine{T}}) where {T<:Number}
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
# appending, prepending , deleting and splitting at

"prepend `Point`s to a MLine"
function prepend!(m::MLine{T},p::Vector{Point{T}}) where T
    prepend!(m.v,p)
    reconfigure!(m)
end

function reconfigure!(m::MLine)
    # after having updated some objects, need to recompute n
    m.xvec = [i.x for i in m.v]
    m.n = length(m.v)
    m.xrange = extrema(m.xvec)
    m.yrange = (min_y(m),max_y(m))
end

"delete an index"
function delete!(m::MLine,idx)
    deleteat!(m.v,idx)
    reconfigure!(m)
end

"append points to a MLine"
function append!(m::MLine{T},p::Vector{Point{T}}) where T
    append!(m.v,p)
    reconfigure!(m)
end

"insert a single value at an interior index"
function insert!(m::MLine{T},v::Point{T},idx::Int) where T
    insert!(m.v,idx,v)
    reconfigure!(m)
end

"""
    splitat(m::MLine,j::Int,repeat_boundary::Bool=true)

Splits a `MLine` object after given index and returns 2 new `MLine`s as a tuple. If `repeat_boundary` is true, then the separating index is the first point of the second new `MLine`.
"""
function splitat(m::MLine,j::Int,repeat_boundary::Bool=true)
    m1 = MLine(m.v[1:j])
    if repeat_boundary
        m2 = MLine(m.v[j:end])
    else
        m2 = MLine(m.v[j+1:end])
    end
    return (m1,m2)
end

"sort a `MLine` along x-grid"
function sortx!(m::MLine)
    sort!(m.v)  # sorts v
    reconfigure!(m)
end


"""
    intersect(L1::MLine,L2::MLine,s::Int)

Intersecting two lines returns a [`Point`](@ref)

Both L1 and L2 are lines with identical support.

return (Point,true) where true indicates that this point should be added to the envelope
"""
function intersect(L1::MLine,L2::MLine,s::Int)
    xlo,vlo = (L1[s].x,L1[s].y)
    xhi,vhi = (L2[s+1].x,L2[s+1].y)

    # println("xlo,vlo = $xlo,$vlo")
    # println("xhi,vhi = $xhi,$vhi")

    tf = false

    # check simple cases: boundaries?
    if L1[s] == L2[s]
        return (L1[s],false)
    elseif L1[s+1] == L2[s+1]
        return (L1[s+1],false)
    else
        t1 = interp(L1,[xhi])[1].y
        t2 = interp(L2,[xlo])[1].y
        if !isfinite(t1)
            # println("t1 = $t1")
            x_x = xlo
            v_x = vlo
        elseif !isfinite(t2)
            # println("t2 = $t2")
            x_x = xhi
            v_x = vhi
        else
            f_closure(z) = interp(L1,[z])[1].y - interp(L2,[z])[1].y
            if f_closure(xlo) * f_closure(xhi) > 0
                # not opposite signs, no zero to be found
                x_x = xlo
                v_x = vhi
            else
                x_x = fzero(f_closure,xlo,xhi)
                v_x = interp(L1,[x_x])[1].y
                if v_x == typemin(eltype(L2.v[1])) || v_x == NaN
                    x_x = xlo
                    v_x = vhi
                else
                    tf = true
                end
            end
        end
        return (Point(x_x,v_x),tf)
    end
end

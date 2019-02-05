
"""

## Envelope

Holds an array of `Line`s, the upper envelope of those lines, and a vector of `Point`s marking the intersections between lines.

### Fields

* `L      `: Vector of [`Line`]@ref
* `env    `: The upper envelope (i.e pointwise maximum) over `L`, itself a [`Line`]@ref
* `env_set`: `true` if envelope vector was set.
* `isects `: Vector of intersections between `Line`s in `L`
* `removed`: Vector of Points removed from each Line `env` during assembly

"""
mutable struct Envelope{T<:Number}
    L      :: Vector{Line{T}}
    env    :: Line{T}
    env_set :: Bool
    isects :: Vector{Point{T}}
    removed :: Vector{Vector{Int}}
    vbound :: T
    function Envelope(x::T) where {T<:Number}
        this = new{T}()
        this.L = Line{T}[]
        this.env = Line([typemin(T)],[typemin(T)])
        this.env_set = false
        this.isects = Point{T}[]
        this.removed = Vector{Vector{Int}}[ ]
        this.vbound = zero(T)
        return this
    end
    function Envelope(e::Line{T}) where {T<:Number}
        this = new{T}()
        this.L = Line{T}[]
        this.env = deepcopy(e)
        this.env_set = true
        this.isects = Point{T}[]
        this.removed = Vector{Vector{Int}}[ ]
        this.vbound = zero(T)
        return this
    end
    function Envelope(l::Vector{Line{T}}) where {T<:Number}
        this = new{T}()
        this.L = deepcopy(l)
        this.env = Line([typemin(T)],[typemin(T)])
        this.env_set = false
        this.isects = Point{T}[]
        this.removed = Vector{Vector{Int}}[ ]
        this.vbound = zero(T)
        return this
    end
end
function show(io::IO, ::MIME"text/plain", en::Envelope{T}) where {T<:Number}
    print(io,"$T Envelope\n")
    print(io,"env Line set?: $(en.env_set) \n")
    print(io,"num of `Line`s: $(length(en.L))\n")
    print(io,"num of intersections: $(length(en.isects))\n")
    print(io,"num of pts removed: $(length(en.removed))\n")
end
function show(io::IO, en::Envelope{T}) where {T<:Number}
    print(io,length(en.env),"-point $T Envelope")
end

size(e::Envelope) = size(e.L)
eltype(e::Envelope) = eltype(e.L) 
bound(e::Envelope) = e.vbound
getx(en::Envelope) = en.env.x
gety(en::Envelope) = en.env.y
gets(en::Envelope) = en.isects
getr(en::Envelope) = en.removed
# getLine(en::Envelope,j::Int) = en.L[j]
# function set!(en::Envelope{T},L::Line{T}) where {T<:Number}
#     en.env = L
# end
# function set!(en::Envelope{T},id::Int,l::Line{T}) where {T<:Number}
#     en.L[id] = l
# end
# Base.setindex!(en::Envelope,l::Line,id::Int) = en.L[id] = l
Base.getindex(en::Envelope,id::Int) = en.L[id]



"""
    removed!(e::Envelope)

Find which points from each `Line` did not end up in the 
`env` and write them to `removed`.
"""
function removed!(e::Envelope)
    if !e.env_set
        error("you need to set an `upper_env!` first.")
    end
    for l in 1:length(e.L)
        # version 0.7 uses
        # findall((!in)(b),a)
        nix = map(x->!in(x,getx(e)),e.L[l].x)
        niy = map(x->!in(x,gety(e)),e.L[l].y)
        push!(e.removed, find(nix .| niy))
    end
end

function remove_c!(ve::Envelope,ce::Envelope)
    for il in 1:length(ve.L)
        if length(ve.removed[il]) > 0
            del = findin(ce.L[il].x,[i.x for i in ve.removed[il]])
            delete!(ce.L[iL],del)
        end
    end
end


"""
    splitLine(m::Line)

splits a `Line` object at wrong EGM solution points. Wrong solutions appear in kinked regions.
"""
function splitLine(o::Line{T}) where T<:Number

    # 1) find all jump-backs in x-grid
    # println(o.x)
    ii = o.xvec[2:end].>o.xvec[1:end-1]  
    # info("splitLine: ii = $(find(.!(ii)))")
    # info("splitLine: x = $(o.x[find(.!(ii))])")

    # 2) if no backjumps at all, exit
    if all(ii)  
        # return as an Envelope
        return Envelope(o)
    else
    # 3) else, identify subsets
        i = 1
        sections = Line{T}[]  # an array of Lines
        new_sections = Line{T}[]  # an array of Lines
        while true
            # println(ii)
            j = findfirst(ii .!= ii[1])  # identifies all indices within kinked region from left to right until the first kink
            # println(j)

            # if no more kinks
            if j==0
                if i > 1
                    # add remaining Line
                    push!(sections,o)
                end
                # then break
                break
            end
            newm,o = splitat(o,j)  # split old Line at j 
            push!(sections,newm)
            ii = ii[j:end] # chop off from jump index
            i += 1
        end

        # all the ones with 2 un-sorted x corrdinates are illegal lines from connection of two proper ones
        # discard those
        # l2 = [length(s.x)==2 for s in sections]
        ns = [!issorted(s.x) && length(s.x)==2 for s in sections]

        # 4) get rid of illegal sections on x
        # e = Envelope(0.0)
        # for s in eachindex(sections)
        #     if issorted(sections[s])
        #         push!(e.L,sections[s])
        #     else
        #         push!(e.removed,[Point(l.x[jx],l.y[jx],i=jx) for jx in ix_nenv])
        #         push!(blacklist,sections[s].x)
        #     end
        # end

        for s in eachindex(sections)
            if !ns[s]
                if !issorted(sections[s])
                    sort!(sections[s])
                end
                push!(new_sections,sections[s])
            end
        end
        return Envelope(new_sections)
    end
end

function upper_env!(e::Envelope{T}) where T<:Number
    # 5) compute upper envelope of all lines
        # - get all x's from all s and sort into a vector xx
        # - interpolate(extrapolate) all s on xx
        # - how to deal with points at which some Line is infeasible?

    if length(e.L)<2
        error("an upper envelope requires by definition at least 2 lines.")
    end

    # - get all x's from all Lines and sort into a vector xx
    xx = sort(unique(vcat([l.xvec for l in e.L]...)))
    n = length(xx)

    # - interpolate(extrapolate) all Ls on xx
    # this returns a matrix (length(L),n)
    # i.e. each row is the interpolation 
    yy = interp(e.L,xx)

    # find the top line at each point in xx
    r_idx = linemax(yy)  # get row indices only: the row index tells us which Line was optimal at that point.

    env = Line([yy[r_idx[i]][i] for i in 1:length(r_idx)]) # envelope over all lines: just pick max points from each line

    # Identify changes in optimal Line
    # switch in top line after index s (indexing global support xx)
    # s tells us after which position in xx we have a change in optimal line
    s = find(r_idx[2:end].!=r_idx[1:end-1])
    isec = Point[]
    isec_s = Int[]

    if length(s) > 0
        println("r_idx = $r_idx")
        println("s = $s")
        # add intersection points between lines
        # an intersection occurs after index s
        for i in s
            if isapprox(yy[r_idx[i]][i],yy[r_idx[i+1]][i])
                # same point:
                # push!(isec,yy[r_idx[i]][i])
            elseif isapprox(yy[r_idx[i]][i+1],yy[r_idx[i+1]][i+1])
                # same point:
                # push!(isec,yy[r_idx[i]][i+1])
            else
                push!(isec_s, i) # record index position
                push!(isec,intersect(yy[r_idx[i]],yy[r_idx[i+1]], i ))
            end
        end
        # add isecs to env
        if length(isec) > 0
            ic = 0
            for i in 1:length(isec)
                insert!(env,isec[i],isec_s[i] + ic)
                ic += 1  # if you already inserted something before, need to shift indices
            end
        end
    end

    e.env = env 
    e.isects = isec
    @assert(issorted(e.env))
    # say that you have now set an upper envelope on this object
    e.env_set = true
    return nothing
end


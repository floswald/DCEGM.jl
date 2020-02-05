
"""

## Envelope

Holds an array of `MLine`s, the upper envelope of those MLines, and a vector of `Point`s marking the intersections between MLines.

### Fields

* `L      `: Vector of [`MLine`]@ref
* `env    `: The upper envelope (i.e pointwise maximum) over `L`, itself a [`MLine`]@ref
* `env_set`: `true` if envelope vector was set.
* `isects `: Vector of intersections between `MLine`s in `L`
* `removed`: Vector of indices of Points removed from each MLine `env` during assembly
* `vbound` : Value on lower bound of asset domain

"""
mutable struct Envelope{T<:Number}
    L      :: Vector{MLine{T}}
    env    :: MLine{T}
    env_set :: Bool
    isects :: Vector{Point{T}}
    removed :: Vector{Vector{Int}}
    vbound :: T
    # Envelope(1) builds a test object with env_set = false
    function Envelope(x::T) where {T<:Number}
        this = new{T}()
        this.L = MLine{T}[]
        this.env = MLine([typemin(T)],[typemin(T)])
        this.env_set = false
        this.isects = Point{T}[]
        this.removed = Vector{Vector{Int}}[ ]
        this.vbound = zero(T)
        return this
    end
    # Envelope(MLine) builds an object with the MLine as THE envelope
    function Envelope(e::MLine{T}) where {T<:Number}
        this = new{T}()
        this.L = MLine{T}[]
        this.env = deepcopy(e)
        this.env_set = true
        this.isects = Point{T}[]
        this.removed = Vector{Vector{Int}}[ ]
        this.vbound = zero(T)
        return this
    end
    # Envelope(Vector{MLine{T}}) are constituting lines, but no env_set yet
    function Envelope(l::Vector{MLine{T}}) where {T<:Number}
        # println("calling Vector{MLine{T}} constr now")
        # println(l[1].v)
        # println(l[2].v)
        this = new{T}()
        this.L = deepcopy(l)
        this.env = MLine([typemin(T)],[typemin(T)])
        this.env_set = false
        this.isects = Point{T}[]
        this.removed = Vector{Vector{Int}}[ ]
        this.vbound = zero(T)
        return this
    end
end
function show(io::IO, ::MIME"text/plain", en::Envelope{T}) where {T<:Number}
    print(io,"$T Envelope\n")
    print(io,"env MLine set?: $(en.env_set) \n")
    print(io,"num of `MLine`s: $(length(en.L))\n")
    print(io,"num of intersections: $(length(en.isects))\n")
    print(io,"num of pts removed: $(length(en.removed))\n")
end
function show(io::IO, en::Envelope{T}) where {T<:Number}
    print(io,length(en.env),"-point $T Envelope")
end

size(e::Envelope) = size(e.L)
eltype(e::Envelope) = eltype(e.L)
bound(e::Envelope) = e.vbound
# getx(en::Envelope) = en.env.x
# gety(en::Envelope) = en.env.y
gets(en::Envelope) = en.isects
getr(en::Envelope) = en.removed
# getMLine(en::Envelope,j::Int) = en.L[j]
# function set!(en::Envelope{T},L::MLine{T}) where {T<:Number}
#     en.env = L
# end
# function set!(en::Envelope{T},id::Int,l::MLine{T}) where {T<:Number}
#     en.L[id] = l
# end
# Base.setindex!(en::Envelope,l::MLine,id::Int) = en.L[id] = l
Base.getindex(en::Envelope,id::Int) = en.L[id]



"""
    removed!(e::Envelope)

Find which points from each `MLine` did not end up in the
`env` and write them to `removed`.
"""
function removed!(e::Envelope)
    if !e.env_set
        error("you need to set an `upper_env!` first.")
    end
    for l in 1:length(e.L)
        push!(e.removed,findall(.!(e.L[l].v .∈ Ref(e.env.v))))   # find indices in MLine L[l] that are not in
    end
end


"""
    to_remove_c(ve::Envelope,ce::Envelope)

Given a value function [`Envelope`](@ref) which has a list of indices of `removed` points,
this function returns indices to remove from the policy function
"""
function to_remove_c(ve::Envelope,ce::Envelope)
    cx = getx(ce.env.v)  # x values in consumption function
    rmidx = Vector{Int}[]

    for il in 1:length(ve.L)
        if length(ve.removed[il]) > 0
            rmv = ve.L[il].v[ve.removed[il]]  # (x,y) values removed from env
            push!(rmidx, findall(cx .∈ Ref(getx(rmv))))  # consumption indices to remove
        end
    end
    # if length(rmidx) > 0
    #     deleteat!(ce.env.v, unique(sort(vcat(rmidx...))))
    # end
    unique(sort(vcat(rmidx...)))
end


"""
    splitLine(m::MLine)

splits a `MLine` object at wrong EGM solution points. Wrong solutions appear in kinked regions.
"""
function splitLine(o::MLine{T}) where T<:Number

    # 1) find all jump-backs in x-grid
    # println(o.x)
    xvec = getx(o)
    ii = xvec[2:end].>xvec[1:end-1]
    # info("splitLine: ii = $(find(.!(ii)))")
    # info("splitLine: x = $(o.x[find(.!(ii))])")

    # 2) if no backjumps at all, exit
    if all(ii)
        # return as an Envelope
        return Envelope(o)
    else
    # 3) else, identify subsets
        i = 1
        sections = MLine{T}[]  # an array of Lines
        new_sections = MLine{T}[]  # an array of Lines
        while true
            # println(ii)
            j = findfirst(ii .!= ii[1])  # identifies all indices within kinked region from left to right until the first kink
            # println(j)

            # if no more kinks
            if isnothing(j)
                if i > 1
                    # add remaining MLine
                    push!(sections,o)
                end
                # then break
                break
            end
            newm,o = splitat(o,j)  # split old MLine at j 
            push!(sections,newm)
            ii = ii[j:end] # chop off from jump index
            i += 1
        end

        # all the ones with 2 un-sorted x corrdinates are illegal lines from connection of two proper ones
        # discard those
        # l2 = [length(s.v)==2 for s in sections]
        # ns = [!issorted(s.v) && length(s.v)==2 for s in sections]

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

        for s in sections
            if (!issorted(s.v)) && (length(s.v)==2)
                # disregard those
            else
                if !issorted(s)
                    sortx!(s)
                end
                push!(new_sections,s)
            end
        end
        if length(new_sections) == 1
            println("gotcha")
            println(new_sections[1].v)
            println("sections[1].v = $(sections[1].v)")
            println("sections[2].v = $(sections[2].v)")
            # plot(Envelope(sections))

            # savefig("check.png")
            println("xvec = $xvec")
            println("ii = $ii")
            error()
            gui()
            # error()
        end
        return Envelope(new_sections)
    end
end


"""
    upper_env!(e::Envelope{T}; do_intersect::Bool=false) where T<:Number

### Outline

This function computes the *upper envelope* over it's constituting lines.
In particular:

1. generates a common support `xx` by concatenating the `x` coords of each `MLine` in `e.L`
2. interpolates all `MLine` in `e.L` on that support `xx`
3. finds the index of the `MLine` with maximal `y`-value for each `xx`
4. assembles the *upper envelope* by just combining all values from 3. into a new `MLine`.

### Optional: `do_intersect`

By setting keyword arg `do_intersect = true`, the function will identify the indices in `xx` after which a change of optimal line occurs. Then it will proceed to find the precise *intersection* between the two constituting `MLine`s involved in this change.

"""
function upper_env!(e::Envelope{T}; do_intersect::Bool=false) where T<:Number
    # 5) compute upper envelope of all lines
        # - get all x's from all s and sort into a vector xx
        # - interpolate(extrapolate) all s on xx
        # - how to deal with points at which some MLine is infeasible? => set to -Inf

    if length(e.L)<2
        println(e.env.v)
        println(e.L[1].v)
        println(e.removed)
        error("an upper envelope requires by definition at least 2 lines.")
    end

    # - get all x's from all Lines and sort into a vector xx
    xx = sort(unique(reduce(vcat,[getx(l) for l in e.L])))
    n = length(xx)

    # - interpolate(extrapolate) all Ls on xx
    # this returns an array (length(L)) of interpolated MLines
    yy = interp(e.L,xx)  # by default no extrapolation here: if out of domain, MLine is -Inf
    # println("yy[1] = $(yy[1].v)")
    # println("yy[2] = $(yy[2].v)")
    # println("yy[3] = $(yy[3].v)")
    # find the top line at each point in xx
    r_idx = linemax(yy)  # get row indices only: the row index tells us which MLine was optimal at that point.

    # envelope over all lines: just pick max points for each xx
    env = MLine([yy[r_idx[i]][i] for i in 1:length(r_idx)])


    if do_intersect
        # Identify changes in optimal MLine
        # switch in top line after index s (indexing global support xx)
        # s tells us after which position in xx we have a change in optimal line
        s = findall(r_idx[2:end].!=r_idx[1:end-1])

        isec = Point[]
        isec_s = Int[]

        if length(s) > 0
            # println("r_idx = $r_idx")
            # println("s = $s")
            # add intersection points between lines
            # an intersection occurs after index s
            ioff = 0    # offsetting add indices
            joff = 1    # isec indices
            for i in s
                tmp = intersect(yy[r_idx[i]],yy[r_idx[i+1]], i )
                if !isnothing(tmp)
                    push!(isec_s, i) # record index position
                    push!(isec,tmp[1])
                    if tmp[2]
                        insert!(env,isec[joff],isec_s[joff] + ioff)
                        ioff += 1  # added additional index: need to shift right now!
                    end
                    joff += 1
                end
            end
        end

        e.isects = isec
    end
    e.env = env
    e.env_set = true
    return nothing
end

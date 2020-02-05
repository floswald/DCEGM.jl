
"""

## Envelope

Holds an array of `MLine`s, the upper envelope of those MLines, and a vector of `Point`s marking the intersections between MLines.

### Fields

* `L      `: Vector of [`MLine`]@ref
* `dirty  `: Initial, potentially backwards bending upper envelope
* `env    `: Cleaned upper envelope
* `env_clean`: `true` if clean envelope vector was set.
* `isects `: Vector of intersecting `Point`s (new points in clean env)
* `removed`: Vector of indices of Points removed from original MLine during assembly
* `vbound` : Value on lower bound of asset domain

"""
mutable struct Envelope{T<:Number}
    L         :: Vector{MLine{T}}
    env       :: MLine{T}
    dirty     :: MLine{T}
    env_clean :: Bool
    isects    :: Vector{Point{T}}
    removed   :: Vector{Int}
    vbound    :: T
    # Envelope(1) builds a test object 
    function Envelope(x::T) where {T<:Number}
        this = new{T}()
        this.L = MLine{T}[]
        this.dirty = MLine([typemin(T)],[typemin(T)])
        this.env = MLine([typemin(T)],[typemin(T)])
        this.env_clean = false
        this.isects = Point{T}[]
        this.removed = Int[ ]
        this.vbound = zero(T)
        return this
    end
    # Envelope(MLine) builds an object with the MLine as dirty envelope
    function Envelope(d::MLine{T}) where {T<:Number}
        this = new{T}()
        this.L = MLine{T}[]
        this.dirty = deepcopy(d)
        this.env = MLine([typemin(T)],[typemin(T)])
        this.env_clean = false
        this.isects = Point{T}[]
        this.removed = Int[ ]
        this.vbound = zero(T)
        return this
    end
    function Envelope(d::MLine{T},c::MLine{T}) where {T<:Number}
        this = new{T}()
        this.L = MLine{T}[]
        this.dirty = deepcopy(d)
        this.env = deepcopy(c)
        this.env_clean = true
        this.isects = Point{T}[]
        this.removed = Int[ ]
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
        this.dirty = MLine([typemin(T)],[typemin(T)])
        this.env = MLine([typemin(T)],[typemin(T)])
        this.env_clean = false
        this.isects = Point{T}[]
        this.removed = Int[ ]
        this.vbound = zero(T)
        return this
    end
end
function show(io::IO, ::MIME"text/plain", en::Envelope{T}) where {T<:Number}
    print(io,"$T Envelope\n")
    print(io,"cleaned env set?: $(en.env_clean) \n")
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



# """
#     removed!(e::Envelope)

# Find which points from each `MLine` did not end up in the
# `env` and write them to `removed`.
# """
# function removed!(e::Envelope)
#     if !e.env_set
#         error("you need to set an `upper_env!` first.")
#     end
#     for l in 1:length(e.L)
#         push!(e.removed,findall(.!(e.L[l].v .∈ Ref(e.env.v))))   # find indices in MLine L[l] that are not in envelope
#     end
# end

"""
    removed!(e::Envelope,L::MLine)

Find which points from initial `MLine` did not end up in the
`env` and write them to `removed`.
"""
function removed!(e::Envelope,L::MLine)
    if !e.env_set
        error("you need to compute an `upper_env!` first.")
    end
    push!(e.removed,findall(.!(L.v .∈ Ref(e.env.v))))   # find indices in L that are not in final envelope
    return nothing
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

        for s in sections
            # fedor keeps those
            # if (!issorted(s.v)) && (length(s.v)==2)
                # disregard those
            # else
                if !issorted(s)
                    sortx!(s)
                end
                push!(new_sections,s)
            # end
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
3. Remembers which indices of each line have been extrapolated, so they cannot be the optimum
4. finds the index of the `MLine` with maximal `y`-value for each `xx` (among non extrapolated ones)
5. assembles the *upper envelope* by just combining all values from 3. into a new `MLine`. This means to collect intersections between lines (i.e. new points) along the way as well.
"""
function upper_env!(e::Envelope{T}) where T<:Number
    # compute upper envelope of all lines
        # - get all x's from all s and sort into a vector xx
        # - interpolate(extrapolate) all s on xx
        # - how to deal with points at which some MLine is infeasible? => set to -Inf

    if length(e.L)<2
        error("an upper envelope requires by definition at least 2 lines.")
    end

    # - get all x's from all Lines and sort into a vector xx
    xx = sort(unique(reduce(vcat,[getx(l) for l in e.L])))

    # - interpolate(extrapolate) all Ls on xx
    # this returns an array (length(L)) of interpolated MLines
    yy = interp(e.L,xx,extrap=false)  # set a line to -Inf if xx not in its domain

    n = length(xx)   # consider all points

    r_idx = linemax(yy)  # get row indices only: the row index tells us which MLine was optimal at that point.
    # envelope over all lines: just pick max points for each xx
    top = [yy[r_idx[i]][i] for i in 1:length(r_idx)]

    # build up upper envelope one by one point in xx
    # ==============================================

    env = Point{T}[]
    push!(env,top[1])  # first point in envelope point vector
    isec = Point[]     # empty vector of point to collect intersections

    k0 = r_idx[1]  # start index
    for j in 2:n
        k1 = r_idx[j]
        if k0 != k1
            # @debug "switching from $k0 to $k1 at index $j"
            # is there an intersection between both lines?
            ln1 = k0; ln2 = k1  # candidate line indices
            L1 = e.L[ln1]; L2 = e.L[ln2]  # candidate lines
            pt1 = xx[j-1]; pt2 = xx[j]  # intersection ∈ [pt1,pt2]

            y1 = interp(e.L[ln1],[pt1,pt2],extrap=true)
            y2 = interp(e.L[ln2],[pt1,pt2],extrap=true)

            # if intersection is *on* either pt1 or pt2, those are intersections
            onboth = y2.v .== y1.v
            # debug> y2.v .== y1.v
            # 2-element BitArray{1}:
            #  1
            #  0
            if any(onboth)
                @assert sum(onboth) == 1
                push!(env,y2.v[onboth][1])
                push!(isec,y2.v[onboth][1])
            else
                # check neither y1 nor y2 are extrapolated in both points
                notboth = (y1.iextrap != [1,2]) && (y2.iextrap != [1,2])
                # and that we have different signs on both ends of search interval
                f_closure(z) = interp(L1,[z])[1].y - interp(L2,[z])[1].y
                diffsign = f_closure(pt1) * f_closure(pt2) < 0

                if notboth && diffsign

                    # find intersection point
                    while true
                        # @debug "checking lines ln1=$ln1 ln2=$ln2 in [$(round(pt1,digits=3)),$(round(pt1,digits=3))]"
                        pt3 = fzero(x -> interp(e.L[ln1],[x])[1].y - interp(e.L[ln2],[x])[1].y,pt1,pt2)
                        y3  = interp(e.L[ln1],[pt3])[1].y  # get function value of ln1 at intersection
                        # @debug "found intersection" pt3=pt3 y3=y3

                        # are there other lines above this intersection point?
                        # interpolate *all* lines in new point pt3
                        yy2 = interp(e.L,[pt3],extrap=false)
                        y3, ln3 = findmax([yy2[i].v[1].y for i in 1:length(e.L)])

                        if ln3 == ln2 || ln3 == ln1
                            # @debug "no higher line at" pt3=pt3
                            # no new lines above at pt3!
                            # add intersection point to the end of the env points
                            push!(env, Point(pt3,y3))
                            # add intersection to intersections container
                            push!(isec, Point(pt3,y3))
                            # additional intersections before next point?
                            if ln2 == k1
                                # no other intersection.
                                # if ln2 is new optimal line and there is no ln3 above,
                                # there cannot be an additional intersection to the right of pt3.
                                break

                            else
                                #
                                # @debug "additional iscecs. updating" ln1=ln2 pt1=pt3 ln2=k1 pt2=xx[j]
                                ln1=ln2
                                pt1=pt3
                                ln2=k1
                                pt2=xx[j]

                            end

                        else
                            # there is indeed a third line over point pt3
                            # what we have is not the upper envelope yet.
                            # redo search.
                            # there must be an intersection to the left of pt3
                            # @debug "additional isec to left. updating" ln2=ln3 pt2=pt3
                            ln2 = ln3 # new candidate line's index
                            pt2 = pt3 # new right boundary of search interval

                        end
                    end
                else
                    # @debug "no valid intersection" notboth=notboth diffsign=diffsign
                end
            end  # if candidate intersection on both lines
        end # if k0 != k1
        # regular addition of points
        # if current xx[j] is a member of line k1, add the point
        # if j is the last index, add the point (no intersections to the right)
        if (xx[j] ∈ e.L[k1].v) || j==n
            # @debug "adding regular point" x=xx[j] Lx=getx(e.L[k1].v)
            # look how beautiful that is: x ∈ v where v is an array of x-y points
            push!(env, top[j])  # push that point onto the envelope
        end
        k0 = k1 # update current index
    end # end j in 2:n
    e.env = MLine(env)
    e.env_set = true
    e.isects = isec
    return nothing
end


"""
    secondary_envelope(L::MLine)

* computes secondary envelope of a `MLine`
* splits line at backward bends and computes upper envelope over those
* returns an `Envelope` object with fields `removed` giving the indices of removed points for each constituting line.
"""
function secondary_envelope(L::MLine)

    # identify loop-backs and split line
    e = splitLine(L)   

    if !e.env_set
        # compute upper envelope
        upper_env!(e)
    end

    # record which points in original L did not make it into final e.env
    removed!(L,e)

    return e

end

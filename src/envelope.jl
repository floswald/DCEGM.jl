
"""

## Envelope

Holds an array of `Line`s, the upper envelope of those lines, and a vector of `Point`s marking the intersections between lines.

### Fields

* `L      `: Vector of [`Line`]@ref
* `env    `: The upper envelope (i.e pointwise maximum) over `L`, itself a [`Line`]@ref
* `env_set`: `true` if envelope vector was set.
* `isects `: Vector of intersections between `Line`s in `L`
* `removed`: Vector of Points removed from `env` during assembly

"""
mutable struct Envelope{T<:Number}
    L      :: Vector{Line{T}}
    env    :: Line{T}
    env_set :: Bool
    isects :: Vector{Point{T}}
    removed :: Vector{Vector{Point{T}}}
    vbound :: T
    function Envelope(x::T) where {T<:Number}
        this = new{T}()
        this.L = Line{T}[]
        this.env = Line([typemin(T)],[typemin(T)])
        this.env_set = false
        this.isects = Point{T}[]
        this.removed = Vector{Point{T}}[Point{T}[] ]
        this.vbound = zero(T)
        return this
    end
    function Envelope(e::Line{T}) where {T<:Number}
        this = new{T}()
        this.L = Line{T}[]
        this.env = deepcopy(e)
        this.env_set = true
        this.isects = Point{T}[]
        this.removed = Vector{Point{T}}[Point{T}[] ]
        this.vbound = zero(T)
        return this
    end
    function Envelope(l::Vector{Line{T}}) where {T<:Number}
        this = new{T}()
        this.L = deepcopy(l)
        this.env = Line([typemin(T)],[typemin(T)])
        this.env_set = false
        this.isects = Point{T}[]
        this.removed = Vector{Point{T}}[Point{T}[] ]
        this.vbound = zero(T)
        return this
    end
end
function show(io::IO, ::MIME"text/plain", en::Envelope{T}) where {T<:Number}
    print(io,"$T Envelope\n")
    print(io,"env Line set?: $(en.env_set) \n")
    print(io,"num of `Line`s: $(length(en.L))\n")
    print(io,"num of intersections: $(length(en.isects))\n")
    print(io,"num of pts removed: $(sum(length(i) for i in en.removed))\n")
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
    for l in e.L
        ix = findin(getx(e),l.x) # find indices in e.x of values that appear (somewhere) in l
        iy = findin(gety(e),l.y) # find indices in e.y of values that appear (somewhere) in l
        # if ix == iy , (ix,iy) is a valid point
        points = findin(ix,iy) 
        iix = ix[points]
        # iiy = iy[points]
        nx = setdiff(l.x,getx(e)[iix])   # x points l but not in env
        ny = setdiff(l.y,gety(e)[iix])   # y points l but not in env
        @assert(length(nx)==length(ny))
        if length(nx) > 0
            push!(e.removed,[Point(nx[jx],ny[jx]) for jx in 1:length(nx)])
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
    ii = o.x[2:end].>o.x[1:end-1]  
    info("splitLine: ii = $(find(.!(ii)))")
    info("splitLine: x = $(o.x[find(.!(ii))])")

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
    xx = sort(unique(vcat([l.x for l in e.L]...)))
    n = length(xx)

    # - interpolate(extrapolate) all Ls on xx
    # this returns a matrix (length(L),n)
    # i.e. each row is the interpolation 
    yy = interp(e.L,xx)

    # find the top line at each point in xx
    val,lin_ind = findmax(yy,1)  # colwise max
    subs = map(x->ind2sub(yy,x),lin_ind)  # get subsript indices of colwise maxima
    r_idx = [i[1] for i in subs] # get row indices only: the row index tells us which Line was optimal at that point.

    # from r_idx we could get which points of each line will be included in the envelope?


    # Identify changes in optimal Line
    # switch in top line after index s (indexing global support xx)
    # s tells us after which position in xx we have a change in optimal line
    s = find(r_idx[2:end].!=r_idx[1:end-1])

    info("upper_env: jumps after $s")
    info("upper_env: jumps at $(xx[s])")



    # Assemble Upper Envelope from Line segments
    # ==========================================

    if length(s)==0
        # there is one complete upper envelope already
        # return
        e.env = Line(xx,yy[r_idx[1],:])
        @assert(issorted(getx(e)))
        # e.isects = [Point(NaN,NaN)]
        # e.removed = [[Point(NaN,NaN)]]
    else
        # sort out which line is top at which index of xx and compute intersection points in between switches
        # s = 1: there is a switch in top line after the first index
        # s = i: there is a switch in top line after the i-th index

        # return Line 
        # The envelope starts with the first line that is on top
        # that line is on top until index s[1] in xx, after which the 
        # top line changes.
        env = Line(xx[1:s[1]],yy[r_idx[s[1]],1:s[1]] )
        isec = [Point(NaN,NaN) for i in 1:length(s)]

        for id_s in eachindex(s)   # id_s indexes line segment in resulting envelope: for n intersecting lines, there are n-1 segments

            js = s[id_s]  # value of index: position in xx
            jx  = subs[js][2]   # colindex of element in yy before switch takes place
            jjx = subs[js+1][2] # colindex of element in yy after switch took place

            # switching from Line to Line
            from = r_idx[js]
            to   = r_idx[js+1]
            info("from L number $(r_idx[js])")
            info("to L Number $(r_idx[js+1])")

            # xx coordinates between which the switching occurs
            # remember xx is a vector as long as size(yy,2)
            x_from = xx[jx ]  # only pick col coordinate
            x_to   = xx[jjx]  # only pick col coordinate
            info("x_from = $(xx[jx ])")  # only pick col coordinate
            info("x_to   = $(xx[jjx])")  # only pick col coordinate

            # values of both lines at those corrdinates
            v_from = yy[subs[js]...]
            v_to   = yy[subs[js+1]...]
            info("v_from = $(yy[subs[js]...])")
            info("v_to   = $(yy[subs[js+1]...])")

            # The intersction of both lines ∈ [x_from,x_to]
            # complication: intersection could also be on the egdes of this interval.
            # this is easy to check:
            if isapprox(yy[from,jx] , yy[to,jx])
                # intersection is on lower bound of interval
                isec[id_s] = Point(x_from,v_from,i=id_s)
                # don't add intersection to envelope!
            elseif isapprox(yy[from,jjx] , yy[to,jjx])
                # intersection is on upper bound of interval
                isec[id_s] = Point(x_to,v_to,i=id_s)
                # don't add intersection to envelope!
            else
                # need to to compute intersection
                f_closure(z) = interp(e.L[to],[z])[1] - interp(e.L[from],[z])[1]
                if f_closure(x_from) * f_closure(x_to) > 0
                    # not opposite signs, no zero to be found
                    x_x = e.L[to]
                    v_x = interp(e.L[to],[x_x])[1]
                else
                    x_x = fzero(f_closure,x_from,x_to)
                    v_x = interp(e.L[to],[x_x])[1]
                    # record intersection as new point
                    isec[id_s] = Point(x_x,v_x,i=id_s,newpoint=true)
                end

                # and add to envelope
                append!(env,x_x,v_x)
                # info("added $(isec[id_s]) to envelope")
                # info("x in env? $(in(x_x,env.x))")
            end

            # sometimes we cannot find an intersection point.
            # either you extrapolate all lines, or you dont
            # 

            # add next line segment to envelope
            # index range s[id_s]+1:s[id_s+1] is
            #     from current switch (next index after): s[id_s]+1
            #     to last index before next switch: s[id_s+1]
            last_ind = id_s==length(s) ? n : s[id_s+1]
            append!(env,xx[js+1:last_ind],yy[to,js+1:last_ind])
        end  # eachindex(s)

        e.env = env 
        e.isects = isec
        # info("x = $(getx(e))")
        for l in e.L
            # collect points that were removed from each line
            ix_env = findin(l.x,getx(e))
            iy_env = findin(l.y,gety(e))
            info("ix_env = $ix_env")
            # info("1:length(l.x) = $(ix_env)")
            ix_nenv = setdiff(1:length(l.x),ix_env)
            # y_nenv = setdiff(l.y,gety(e))
            @assert(issorted(getx(e)))
            info("setdiff(l.x,x) = $ix_nenv")
            # ix = map(x->!in(x,getx(e)),l.x) 
            # iy = map(x->!in(x,gety(e)),l.y) 
            # jj = find(ix .| iy)
            if length(ix_nenv) > 0
                push!(e.removed,[Point(l.x[jx],l.y[jx],i=jx) for jx in ix_nenv])
            end
        end
        @assert(issorted(getx(e)))
        # say that you have now set an upper envelope on this object
        e.env_set = true
        return nothing
    end
end


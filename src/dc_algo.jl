


function minimal_EGM()
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))
    avec          = scaleGrid(0.0,p.a_high,p.na,logorder = 1)
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,title="0 borrowing")
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
        c[it] = vcat(0.0, 1.0 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector. (p.na+1,1)
        m[it] = vcat(0.0, avec .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        plot!(pl,m[it],c[it],label="$it")
    end
    return (m,c,pl)
end

function minimal_EGM_neg(;alow = -1.0,blim=true,brate = false)
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))

    η = NBL(yvec[1],p)  # compute natural borrowing limit: maximal borrowing if c=0 in all future periods and worst income
    if blim
        avec          = [scaleGrid(η[it],p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
        push!(avec, scaleGrid(0.0,p.a_high,p.na,logorder = 1))  # last period
    else
        avec          = [scaleGrid(alow,p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
        push!(avec, scaleGrid(0.0,p.a_high,p.na,logorder = 1))  # last period
    end
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    # m[p.nT]       = [p.a_low,0.0,p.a_high]    # no debt in last period possible
    # c[p.nT]       = [p.cfloor,p.cfloor,p.a_high]
    ti = blim ? "Ntl. br:$brate" : "Fixed. br:$brate"
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,title = ti)
    # cycle back in time
    for it in p.nT-1:-1:1
        if !brate
            w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec[it+1]'   # w1 = y + yshock + R*savings:  next period wealth at all states. (p.ny,p.na)
        else

            w1 = 0.0 .+ exp.(yvec) .+ vcat(avec[it+1][avec[it+1] .< 0 ].*(p.R*1.3),avec[it+1][avec[it+1] .>= 0 ].*(p.R))'  # w1 = y + yshock + R*savings:  next period wealth at all states. (p.ny,p.na)
        end
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        # if it == p.nT-3
            # println("mt+1 = $(m[it+1])")
            # println("wt+1 = $w1")
            # println("ct+1 = $c1")
        # end
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
        c[it] = vcat(0.0, 1 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector ∈ [0,∞]
        m[it] = vcat(avec[it][1], avec[it] .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        # plot!(pl,m[it],c[it],label="$it")
        plot!(pl,m[it],c[it])
    end
    return (m,c,pl)
end

function plot_minimal()
    p1 = minimal_EGM()
    p2 = minimal_EGM_neg()
    p3 = minimal_EGM_neg(brate = true)
    p4 = minimal_EGM_neg(blim = false)
    p5 = minimal_EGM_neg(blim = false,brate = true)
    plot(p1[3],p2[3],p3[3],p4[3],p5[3],layout = (1,5),xlims = (-10,10),size = (800,400))
end


function dc_EGM!(m::FModel,p::Param)

    for it in p.nT:-1:1
        println(it)
        # @info("period = $it")

        if it==p.nT
            for id in 1:p.nD   # work of dont work
                # final period: consume everyting.
                m.c[id,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,p.a_high)) )
                # initialize value function with vf(1) = 0
                m.v[id,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,NaN)) )
                # note that 0.0 as first value of the vfun is not innocuous here!
            end
        else
            for id in 1:p.nD   # current period dchoice
                working = id==1  # working today is id=1
                # what's next period's cash on hand given you work/not TODAY?
                mm1 = m.m1[it][id]
                clamp!(mm1, p.cfloor, Inf)  # dont actullay want to do that: want to allow for neg assets

                # next period consumption and values y-coords
                # for each d-choice
                cmat = fill(-Inf,p.nD,p.na*p.ny)
                vmat = fill(-Inf,p.nD,p.na*p.ny)

                # get future cons and values
                for iid in 1:p.nD
                    c1 = interp(m.c[iid,it+1].env, mm1[:])
                    floory!(c1,p.cfloor)   # floor negative consumption
                    cmat[iid,:] = gety(c1)  # get y-values
                    vmat[iid,:] = vfun(iid,it+1,cmat[iid,:],mm1[:],m.v[iid,it+1],p)
                end

                # get ccp to be a worker
                pwork = working ? ccp(vmat,p) : zeros(size(vmat)[2])

                # get expected marginal utility of that consumption
                mu1 = reshape(pwork .* up(cmat[1,:],p) .+ (1.0 .- pwork) .* up(cmat[2,:],p),p.ny,p.na)

                # get expected marginal value of saving: RHS of euler equation
                # beta * R * E[ u'(c_{t+1}) ]
                RHS = p.beta * p.R *  m.ywgt' * mu1

                # optimal consumption today: invert the RHS of euler equation
                c0 = iup(Array(RHS)[:],p)

                # set optimal consumption function today. endo grid m and cons c0
                cline = MLine(m.avec .+ c0, c0)
                # store
                m.c[id,it] = Envelope(cline)

                # consumption function done.


                # compute value function
                # ----------------------
                if working
                    # ev = reshape(logsum(vmat,p),p.na,p.ny) * m.ywgt[:,iy]
                    ev =  m.ywgt' * reshape(logsum(vmat,p),p.ny,p.na)
                else
                    ev =  m.ywgt' * reshape(vmat[2,:],p.ny,p.na)
                end
                vline = MLine(m.avec .+ c0, u(c0,id==1,p) .+ p.beta * ev[:])
                # vline is an *uncleaned* value function which may have backward bends

                if any(isnan.(ev))
                    println("ev = ")
                    display(ev)
                end



                # vline and cline may have backward-bending regions: let's prune those
                # SECONDARY ENVELOPE COMPUTATION

                if id==1   # only for workers
                    minx = min_x(vline)
                    if minx < vline.v[1].x
                        # non-convex region lies inside credit constraint.
                        # endogenous x grid bends back before the first x grid point.
                        x0 = collect(range(minx,stop = vline.v[1].x,length = max(10,floor(Integer,p.na/10)))) # some points to the left of first x point
                        x0 = x0[1:end-1]
                        y0 = u(x0,working,p) .+ p.beta .* ev[1]
                        prepend!(vline,convert(Point,x0,y0))
                        prepend!(cline,convert(Point,x0,y0))  # cons policy in credit constrained is 45 degree line
                    end

                    # split the vline at potential backward-bending points
                    # and save as Envelope object
                    m.v[id,it] = splitLine(vline)  # splits line at backward bends

                    # if there is just one line (i.e. nothing was split in preceding step)
                    # then this IS a valid envelope
                    # else, need to compute the upper envelope.
                    if !m.v[id,it].env_set
                        upper_env!(m.v[id,it],do_intersect = true)   # compute upper envelope of this
                        # println(m.c[id,iy,it].env)

                        removed!(m.v[id,it])
                        # sortx!(m.c[id,it].env)

                        @assert(issorted(getx(m.v[id,it].env)))
                        # @assert(issorted(getx(m.c[id,it].env)))


                        # if any points were removed from v, need to check c
                        if length(length(getr(m.v[id,it]))) > 0

                            rmidx = to_remove_c(m.v[id,it],m.c[id,it])  # indices to be removed from current c
                        # if length(isecs) > 0
                            insert_left = Point[]
                            insert_right = Point[]

                            # insert new intersections into consumption function
                            isecs = gets(m.v[id,it])
                            consx = getx(m.c[id,it].env)


                            for isec in 1:length(isecs)
                                I = isecs[isec]

                                # interpolate from left
                                jl = findall(consx .< I.x)
                                jl = jl[.!(jl .∈ Ref(rmidx))]  # keep those who are not to be deleted
                                if length(jl) > 0
                                    jl = maximum(jl)  # biggest of those
                                    newleft = MLine(m.c[id,it].env.v[jl:jl+1])
                                    sortx!(newleft)
                                    push!(insert_left, getv(interp(newleft, [ I.x ] ))[1] )
                                else
                                    push!(insert_left,I)
                                end

                                # interpolate from right
                                jr = findall(consx .> I.x)
                                jr = jr[.!(jr .∈ Ref(rmidx))]  # keep those who are not to be deleted
                                if length(jr) > 0
                                    jr = minimum(jr)   # smallest of those
                                    # push!(insert_right, interp(m.c[id,it].env[jr-1:jr], [ I.x ] ) )
                                    newright = MLine(m.c[id,it].env.v[jr-1:jr])
                                    sortx!(newright)
                                    push!(insert_right, getv(interp(newright, [ I.x ] ))[1] )
                                else
                                    push!(insert_right,I)
                                end

                                    # insert two points at I.x into consumption function
                                    # from left, slightly offset by eps() to preserve ordering of x
                                    # from right, exactly on I.x
                                # end
                            end

                            # remove illegal points from c
                            rmidx = to_remove_c(m.v[id,it],m.c[id,it])  # indices to be removed from current c
                            deleteat!(m.c[id,it].env.v, rmidx)
                            # consx = getx(m.c[id,it].env)
                            # unique!(m.c[id,it].env)


                            consx = getx(m.c[id,it].env)
                            # insert new points at intersections
                            for isec in 1:length(insert_left)
                                I = isecs[isec]
                                jr = findfirst(consx .> I.x)
                                if  !isnothing(jr) && jr > 1
                                    insert!(m.c[id,it].env, Point(I.x-5*eps(),insert_left[isec].y) , jr-1)
                                    insert!(m.c[id,it].env, Point(I.x      ,insert_right[isec].y), jr)
                                end
                            end



                            # @assert(issorted(consx))
                            sortx!(m.c[id,it].env)


                            @assert(issorted(getx(m.c[id,it].env)))

                            # now insert all. need second loop because
                            # for isec in 1:length(isecs)
                            #     I = isecs[isec]
                            #     consx = getx(m.c[id,it].env)

                            #     # interpolate from left
                            #     jl = findlast(consx .< I.x)
                            #     push!(insert_left, interp(m.c[id,it].env[jl:jl+1], [ I.x ] ) )

                            #     # interpolate from right
                            #     jr = findfirst(consx .> I.x)
                            #     push!(insert_right, interp(m.c[id,it].env[jr-1:jr], [ I.x ] ) )
                            # end



                            # interpolate from right


                            # add new points twice to accurately describe discontinuity



                                # I = isecs[isec]

                                # if that intersection is a new point
                                # # i.e. intersection was not a member of any `MLine`
                                # if I.new_point
                                #     # insert intersection into env over cons function
                                #     println("I.x = $(I.x)")
                                #     println("I.i = $(I.i)")

                                #     if !issorted(m.c[id,it].env.x)
                                #         println("m.c[id,it].env.x = $(m.c[id,it].env.x)")
                                #         println("m.c[id,it].env.y = $(m.c[id,it].env.y)")
                                #     end

                                #     insert!(m.c[id,it].env,I.x,interp(m.c[id,it].env,[I.x]),I.i)

                                #     # add to both adjacent `MLine` segments:
                                #     # 1) append to end of segment preceding intersection:
                                #     newy = interp(m.c[id,it].L[I.i],[I.x])
                                #     append!(m.c[id,it].L[I.i],I.x,newy)
                                #     # 1) prepend to beginning of segment following intersection:
                                #     prepend!(m.c[id,it].L[I.i+1],I.x,newy)
                                # end
                            # end
                        end
                    end
                else   # if id==1
                    m.v[id,it] = Envelope(vline)
                end

                # store the expected value at the lower boundary
                # in a separate object
                # NOT as the first value in the vfun as Fedor.
                m.v[id,it].vbound = ev[1]

                # this creates the credit constrained region
                prepend!(m.c[id,it].env,[Point(m.avec[1],0.0)])
                prepend!(m.v[id,it].env,[Point(m.avec[1],ev[1])])
                # sortx!(m.c[id,it].env)
                # sortx!(m.v[id,it].env)
                # prepend!(m.v[id,it].env,p.a_low,ev[1])
                # do NOT prepend the value function with the special value from above.
            end # current discrete choice
        end    # if final perio
    end     # loop over time


end





"""
    vfun(id::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Array{Envelope},p::Param)

Calculate the period `it`, discrete choice `id`-specific value function. Avoids interpolation in credit constrained region by using the analytic form of the value function (no need to interpolate expected value function when on the lower bound of assets.)
"""
function vfun(id::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},v::Envelope,p::Param)

    # L = en.L[id]
    # v = en[id,iy,it]

    # computes v_{it}(m) = u(c) + beta v_{it+1}(m1)

    if length(getx(v.env)) < 2
        error("need more than 2 points in envelope object")
    end

    r = fill(NaN,size(m1))
    mask = m1.<getx(v.env)[2]
    mask = it==p.nT ? trues(size(mask)) : mask

    if all(mask)
        # in the credit constrained region:
        r[:] = u(c1,id==1,p) .+ p.beta * bound(v)
    elseif any(mask)
        r[mask] = u(c1[mask],id==1,p) .+ p.beta * bound(v)
        # elsewhere
        r[.!mask] = gety(interp(v.env,m1[.!mask]))
    else
        r[:] = gety(interp(v.env,m1))
    end

    return r
end

"""
    ccp(x::Matrix,p::Param)

Conditional Choice probability of working
"""
function ccp(x::Matrix,p::Param)
    #choice probability of the first row in 2-row matrix
    mx = maximum(x,dims = 1)
    mxx = x.-repeat(mx,size(x)[1],1)   # center values at max for numerical stability
    vec(exp.(mxx[1,:]./p.lambda)' ./ sum(exp.(mxx./p.lambda),dims = 1))
end

"""
    logsum(x::Matrix,p::Param)

Logsum of conditional values used in Expected value function.
"""
function logsum(x::Matrix,p::Param)
    mx = maximum(x,dims = 1)
    mxx = x .- repeat(mx,size(x)[1],1)
    mx .+ p.lambda * log.( sum(exp.(mxx./p.lambda), dims = 1) )
end


"""
    dcegm!(m::GModel,p::Param)

Main body of the DC-EGM algorithm version
"""
function dc_EGM!(m::GModel,p::Param)
    for it in p.nT:-1:1
        println(it)
        # @info("period = $it")

        if it==p.nT
            for iy in 1:p.ny
                for id in 1:p.nD   # work of dont work
                    # final period: consume everyting.
                    # set the consumption function
                    # remember this is a `MLine`, i.e. it has an x and a y Vector
                    # x: endogenous grid m
                    # y: optimal consumption at that grid x
                    m.c[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,p.a_high)) )

                    # initialize value function with vf(1) = 0
                    m.v[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,NaN)) )
                    # note that 0.0 as first value of the vfun is not innocuous here!
                end
            end

        else
            for iy in 1:p.ny
                # @info("iy: $iy")
                for id in 1:p.nD   # current period dchoice
                    working = id==1  # working today is id=1
                    # @info("current period status id: $id")

                    # next period consumption and value y-coords
                    # for each d-choice
                    cmat = fill(-Inf,p.nD,p.na*p.ny)
                    vmat = fill(-Inf,p.nD,p.na*p.ny)

                    if working

                        for iid in 1:p.nD   # next periods' discrete choice!
                            # precomputed next period's cash on hand on all income states
                            # what's next period's cash on hand given you work/not tomorrow?
                            mm1 = m.m1[it+1][iid]

                            # interpolate the iid-choice next period's consumtion function on next cash on hand, given that discrete choice
                            # println(m.c[iid,iy,it+1].env.v)
                            c1 = interp(m.c[iid,iy,it+1].env, mm1[:])
                            floory!(c1,p.cfloor)   # floor negative consumption
                            cmat[iid,:] = gety(c1)  # get y-values

                            vmat[iid,:] = vfun(iid,it+1,cmat[iid,:],mm1[:],m.v[iid,iy,it+1],p)
                        end

                    else  # retirees have no discrete choice - absorbing state

                        iid = 2   # no work next period
                        mm1 = m.m1[it+1][iid]  # non-work wealth
                        c1 = interp(m.c[iid,iy,it+1].env, mm1[:])
                        floory!(c1,p.cfloor)   # floor negative consumption
                        cmat[iid,:] = gety(c1)

                        vmat[iid,:] = vfun(iid,it+1,cmat[iid,:],mm1[:],m.v[iid,iy,it+1],p)

                    end

                    # get ccp to be a worker
                    pwork = working ? ccp(vmat,p) : zeros(size(vmat)[2])

                    # get marginal utility of that consumption
                    mu1 = reshape(pwork .* up(cmat[1,:],p) .+ (1.0 .- pwork) .* up(cmat[2,:],p),p.ny,p.na)
                    # println("mu1 = ")
                    # display(mu1[1:10,:])

                    # get expected marginal value of saving: RHS of euler equation
                    # beta * R * E[ u'(c_{t+1}) | iy ]
                    # need to integrate out Py here
                    RHS = p.beta * p.R *  m.ywgt[iy,:]' * mu1
                    # RHS = p.beta * p.R * mu1 * vec(m.ywgt)
                    # println("RHS = $(RHS[1:10])")

                    # optimal consumption today: invert the RHS of euler equation
                    c0 = iup(Array(RHS)[:],p)

                    # set optimal consumption function today. endo grid m and cons c0
                    cline = MLine((m.avec[it]) .+ c0, c0)
                    # store
                    m.c[id,iy,it] = Envelope(cline)

                    # consumption function done.


                    # compute value function
                    # ----------------------
                    if working
                        # ev = reshape(logsum(vmat,p),p.na,p.ny) * m.ywgt[:,iy]
                        ev =  m.ywgt[iy,:]' * reshape(logsum(vmat,p),p.ny,p.na)
                    else
                        ev =  m.ywgt[iy,:]' * reshape(vmat[2,:],p.ny,p.na)
                    end
                    vline = MLine((m.avec[it]) .+ c0, u(c0,id==1,p) .+ p.beta * ev[:])

                    # println(vline)

                    if any(isnan.(ev))
                        println("ev = ")
                        display(ev)
                    end



                    # vline and cline may have backward-bending regions: let's prune those
                    # SECONDARY ENVELOPE COMPUTATION

                    if id==1   # only for workers
                        minx = min_x(vline)
                        if minx < vline.v[1].x
                            # non-convex region lies inside credit constraint.
                            # endogenous x grid bends back before the first x grid point.
                            x0 = collect(range(minx,stop = vline.v[1].x,length = max(10,floor(Integer,p.na/10)))) # some points to the left of first x point
                            x0 = x0[1:end-1]
                            y0 = u(x0,working,p) .+ p.beta .* ev[1]
                            prepend!(vline,convert(Point,x0,y0))
                            prepend!(cline,convert(Point,x0,y0))  # cons policy in credit constrained is 45 degree line
                        end

                        # split the vline at potential backward-bending points
                        # and save as Envelope object
                        m.v[id,iy,it] = splitLine(vline)  # splits line at backward bends

                        # if there is just one line (i.e. nothing was split in preceding step)
                        # then this IS a valid envelope
                        # else, need to compute the upper envelope.
                        if !m.v[id,iy,it].env_set
                            upper_env!(m.v[id,iy,it])   # compute upper envelope of this
                            # println(m.c[id,iy,it].env)

                            removed!(m.v[id,iy,it])
                            remove_c!(m.v[id,iy,it],m.c[id,iy,it])
                            sortx!(m.c[id,iy,it].env)

                            @assert(issorted(getx(m.v[id,iy,it].env)))
                            # display(hcat(getx(m.v[id,iy,it]),getx(m.c[id,iy,it])))
                            @assert(issorted(getx(m.c[id,iy,it].env)))
                            # insert new intersections into consumption function
                            isecs = gets(m.v[id,iy,it])
                            if length(isecs) > 0
                                for isec in 1:length(isecs)
                                    I = isecs[isec]

                                    # if that intersection is a new point
                                    # # i.e. intersection was not a member of any `MLine`
                                    # if I.new_point
                                    #     # insert intersection into env over cons function
                                    #     println("I.x = $(I.x)")
                                    #     println("I.i = $(I.i)")

                                    #     if !issorted(m.c[id,iy,it].env.x)
                                    #         println("m.c[id,iy,it].env.x = $(m.c[id,iy,it].env.x)")
                                    #         println("m.c[id,iy,it].env.y = $(m.c[id,iy,it].env.y)")
                                    #     end

                                    #     insert!(m.c[id,iy,it].env,I.x,interp(m.c[id,iy,it].env,[I.x]),I.i)

                                    #     # add to both adjacent `MLine` segments:
                                    #     # 1) append to end of segment preceding intersection:
                                    #     newy = interp(m.c[id,iy,it].L[I.i],[I.x])
                                    #     append!(m.c[id,iy,it].L[I.i],I.x,newy)
                                    #     # 1) prepend to beginning of segment following intersection:
                                    #     prepend!(m.c[id,iy,it].L[I.i+1],I.x,newy)
                                    # end
                                end
                            end
                        end
                    else   # if id==1
                        m.v[id,iy,it] = Envelope(vline)
                    end

                    # store the expected value at the lower boundary
                    # in a separate object
                    # NOT as the first value in the vfun as Fedor.
                    m.v[id,iy,it].vbound = ev[1]

                    # this creates the credit constrained region
                    prepend!(m.c[id,iy,it].env,[Point(m.avec[it][1],0.0)])
                    sortx!(m.c[id,iy,it].env)
                    sortx!(m.v[id,iy,it].env)
                    # prepend!(m.v[id,it].env,p.a_low,ev[1])
                    # do NOT prepend the value function with the special value from above.
                end # current discrete choice
            end   # iy
        end    # if final perio
    end     # loop over time
end

function runf()
    p = Param()
    m = FModel(p)
    dc_EGM!(m,p)
    (m,p)
end

function run2()
    p = Param()
    m = GModel(p)
    dc_EGM!(m,p)
    (m,p)
    # plot(m,id=2,xlim=(0,4),ylim=(-100,-20))
    # plot!(m,id=1,linestyle=:dash,label="")
end
function runit()
    p = Param()
    m = GModel(p)
    dc_EGM!(m,p)
    return (m,p)
end
function pp()
    m,p = runit()
    plot(m.v[1,1].env)
end




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

function minimal_EGM_neg(;alow = -1.0,brate = false)
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))

    # η = NBL(yvec[1],p)  # compute natural borrowing limit: maximal borrowing if c=0 in all future periods and worst income
    # if blim
    #     avec          = [scaleGrid(η[it],p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
    #     push!(avec, scaleGrid(0.0,p.a_high,p.na,logorder = 1))  # last period
    # else
        avec          = [scaleGrid(alow,p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
        push!(avec, scaleGrid(0.0,p.a_high,p.na,logorder = 1))  # last period
    # end
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    # m[p.nT]       = [p.a_low,0.0,p.a_high]    # no debt in last period possible
    # c[p.nT]       = [p.cfloor,p.cfloor,p.a_high]
    # ti = blim ? "Ntl. br:$brate" : "Fixed. br:$brate"
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


"""
    dc_EGM!(m::FModel,p::Param)

DCEGM algorithm as in Ishkakov et al.
"""
function dc_EGM!(m::FModel,p::Param)

    for it in p.nT:-1:1
        # println(it)
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
                    ev =  m.ywgt' * reshape(logsum(vmat,p),p.ny,p.na)
                else
                    ev =  m.ywgt' * reshape(vmat[2,:],p.ny,p.na)
                end
                vline = MLine(m.avec .+ c0, u(c0,id==1,p) .+ p.beta * ev[:])

                # SECONDARY ENVELOPE COMPUTATION
                # ==============================

                if id==1   # only for workers
                    minx = min_x(vline)
                    if vline.v[1].x <= minx
                        # normal case - no bending back behind first point
                        # call secondary_envelope on vline
                        m.v[id,it] = secondary_envelope(vline)

                    else
                        # non-convex region lies inside credit constraint.
                        # endogenous x grid bends back before the first x grid point.
                        x0 = collect(range(minx,stop = vline.v[1].x,length = floor(Integer,p.na/10))) # some points to the left of first x point
                        x0 = x0[1:end-1]
                        y0 = u(x0,working,p) .+ p.beta .* ev[1]
                        prepend!(vline,convert(Point,x0,y0))
                        prepend!(cline,convert(Point,x0,x0))  # cons policy in credit constrained is 45 degree line
                        m.c[id,it] = Envelope(cline)
                        m.v[id,it] = secondary_envelope(vline)
                    end

                    # now we have a cleaned value function

                    # if any points were removed from vline:
                    if length(getr(m.v[id,it])) > 0

                        # analyse policy function
                        # =======================

                        rmidx = getr(m.v[id,it])  # indices removed from value function
                        insert_left = Point[]
                        insert_right = Point[]

                        # remove illegal points from c
                        # deleteat!(m.c[id,it].env.v, rmidx)

                        # insert new intersections into consumption function
                        isecs = gets(m.v[id,it])
                        consx = getx(m.c[id,it].env)

                        # compute any new cons-points (intersections in v)
                        # ---------------------------------------

                        for isec in 1:length(isecs)
                            I = isecs[isec]

                            # interpolate from left
                            jl = findall(consx .< I.x)
                            jl = jl[ jl .∉ Ref(rmidx) ]  # keep those who are not to be deleted
                            if length(jl) > 0
                                jl = maximum(jl)  # biggest of those
                                newleft = MLine(m.c[id,it].env.v[jl:jl+1])
                                # @fediskhakov: who guarantees that jl+1 is not to be deleted?
                                # more generally: why do we not delete points before we
                                # insert the new intersections?

                                sortx!(newleft)
                                tmp = getv(interp(newleft, [ I.x ] ))[1]
                                # take this new point at a minimal left shift in x
                                push!(insert_left, Point(tmp.x - 1e3 * eps(), tmp.y) )
                            # else
                            #     push!(insert_left,I)
                            end

                            # interpolate from right
                            jr = findall(consx .> I.x)
                            jr = jr[ jr .∉ Ref(rmidx) ]  # keep those who are not to be deleted
                            if length(jr) > 0
                                jr = minimum(jr)   # smallest of those
                                # push!(insert_right, interp(m.c[id,it].env[jr-1:jr], [ I.x ] ) )
                                newright = MLine(m.c[id,it].env.v[jr-1:jr])
                                sortx!(newright)
                                push!(insert_right, getv(interp(newright, [ I.x ] ))[1] )
                            # else
                            #     push!(insert_right,I)
                            end
                        end # all intersections

                        # remove illegal points from c
                        deleteat!(m.c[id,it].env.v, rmidx)

                        # add new points in twice with a slight offset from left
                        # to preserve the ordering in x.
                        for ileft in 1:length(insert_left)

                            consx = getx(m.c[id,it].env)
                            j = findfirst(consx .> insert_left[ileft].x)  # first point past new intersection
                            insert!(m.c[id,it].env.v,j,insert_left[ileft])  # item is j-th index
                            insert!(m.c[id,it].env.v,j+1,insert_right[ileft])
                        end


                    end
                else   # if id==1
                    m.v[id,it] = Envelope(vline)
                end

                # store the expected value at the lower boundary
                # in a separate object
                m.v[id,it].vbound = ev[1]

                # this creates the credit constrained region
                prepend!(m.c[id,it].env,[Point(m.avec[1],0.0)])
                prepend!(m.v[id,it].env,[Point(m.avec[1],ev[1])])
                # if !issorted(m.c[id,it].env)
                #     println(isecs)
                #     xx = getx(m.c[id,it].env)
                #     ii = findfirst( vcat(0,diff(xx)) .< 0 )
                #     println(m.c[id,it].env.v[ii-1:ii+1])
                    sortx!(m.c[id,it].env)  # sort cons by default
                # end
                # @assert issorted(m.c[id,it].env)
                # sortx!(m.c[id,it].env)
                # sortx!(m.v[id,it].env)
                # prepend!(m.v[id,it].env,p.a_low,ev[1])
                # do NOT prepend the value function with the special value from above.
            end # current discrete choice
        end    # if final perio
    end     # loop over time


end

"""
    dc_EGM!(m::GModel,p::Param)

DCEGM algorithm for a model with state dependence.
"""
function dc_EGM!(m::GModel,p::Param)

    cmat = fill(-Inf,p.nD,p.na)
    vmat = fill(-Inf,p.nD,p.na)
    ctmp = fill(-Inf,p.nD,p.ny,p.na)
    vtmp = fill(-Inf,p.nD,p.ny,p.na)

    for it in p.nT:-1:1
        for iy in 1:p.ny  # current state
            for id in 1:p.nD  # current dchoice
                working = id==1  # working today is id=1

                if it==p.nT
                    # final period: consume everyting.
                    m.c[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,p.a_high)) )
                    # initialize value function with vf(1) = 0
                    m.v[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,NaN)) )
                    # note that 0.0 as first value of the vfun is not innocuous here!
                else
                    # next period consumption and values y-coords
                    # for each d-choice
                    # reset all value matrices
                    fill!(cmat,-Inf)
                    fill!(vmat,-Inf)
                    fill!(ctmp,-Inf)
                    fill!(vtmp,-Inf)

                    for jy in 1:p.ny # future state: owner, renter, income, etc
                        pr = m.ywgt[iy,jy]  # transprob

                        for iid in 1:p.nD  # future dchoice
                            # only feasible choices at this state
                            # if renter, cannot sell

                            m1 = m.m1[it+1][iid][jy,:]  # state specific mvec
                            c1 = interp(m.c[iid,jy,it+1].env, m1) # C(d',y',m')
                            floory!(c1,p.cfloor)   # floor negative consumption
                            ctmp[iid,jy,:] = pr .* gety(c1)
                            vtmp[iid,jy,:] = pr .* vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                        end
                    end # end future state

                    # now get expectations conditional on iy: E[c(t+1,iid,y')|iy]
                    # cmat = integrate over second dimension: nD by na
                    # vmat = integrate over second dimension
                    cmat = dropdims( reduce(+, ctmp, dims = 2), dims = 2)
                    vmat = dropdims( reduce(+, vtmp, dims = 2), dims = 2)

                    # get ccp of choices: P(d'|iy), pwork
                    pwork = working ? ccp(vmat,p) : zeros(size(vmat)[2])

                    # get MU of cons
                    mu1 = pwork .* up(cmat[1,:],p) .+ (1.0 .- pwork) .* up(cmat[2,:],p) # 1,na

                    #RHS
                    RHS = p.beta * p.R * mu1

                    #optimal cons
                    c0 = iup(RHS,p)
                    # set optimal consumption function today. endo grid m and cons c0
                    cline = MLine(m.avec .+ c0, c0)
                    # store
                    m.c[id,iy,it] = Envelope(cline)
                    # consumption function done.

                    # compute value function
                    # ----------------------
                    if working
                        ev =  logsum(vmat,p)
                    else
                        ev =  vmat[2,:]
                    end
                    vline = MLine(m.avec .+ c0, u(c0,id==1,p) .+ p.beta * ev[:])

                    # SECONDARY ENVELOPE COMPUTATION
                    # ==============================

                    if id==1   # only for workers
                        minx = min_x(vline)
                        if vline.v[1].x <= minx
                            # normal case - no bending back behind first point
                            # call secondary_envelope on vline
                            m.v[id,iy,it] = secondary_envelope(vline)

                        else
                            # non-convex region lies inside credit constraint.
                            # endogenous x grid bends back before the first x grid point.
                            x0 = collect(range(minx,stop = vline.v[1].x,length = floor(Integer,p.na/10))) # some points to the left of first x point
                            x0 = x0[1:end-1]
                            y0 = u(x0,working,p) .+ p.beta .* ev[1]
                            prepend!(vline,convert(Point,x0,y0))
                            prepend!(cline,convert(Point,x0,x0))  # cons policy in credit constrained is 45 degree line
                            m.c[id,iy,it] = Envelope(cline)
                            m.v[id,iy,it] = secondary_envelope(vline)
                        end

                        # now we have a cleaned value function

                        # if any points were removed from vline:
                        if length(getr(m.v[id,iy,it])) > 0

                            # analyse policy function
                            # =======================

                            rmidx = getr(m.v[id,iy,it])  # indices removed from value function
                            insert_left = Point[]
                            insert_right = Point[]

                            # remove illegal points from c
                            # deleteat!(m.c[id,iy,it].env.v, rmidx)

                            # insert new intersections into consumption function
                            isecs = gets(m.v[id,iy,it])
                            consx = getx(m.c[id,iy,it].env)

                            # compute any new cons-points (intersections in v)
                            # ---------------------------------------

                            for isec in 1:length(isecs)
                                I = isecs[isec]

                                # interpolate from left
                                jl = findall(consx .< I.x)
                                jl = jl[ jl .∉ Ref(rmidx) ]  # keep those who are not to be deleted
                                if length(jl) > 0
                                    jl = maximum(jl)  # biggest of those
                                    newleft = MLine(m.c[id,iy,it].env.v[jl:jl+1])
                                    # @fediskhakov: who guarantees that jl+1 is not to be deleted?
                                    # more generally: why do we not delete points before we
                                    # insert the new intersections?

                                    sortx!(newleft)
                                    tmp = getv(interp(newleft, [ I.x ] ))[1]
                                    # take this new point at a minimal left shift in x
                                    push!(insert_left, Point(tmp.x - 1e3 * eps(), tmp.y) )
                                # else
                                #     push!(insert_left,I)
                                end

                                # interpolate from right
                                jr = findall(consx .> I.x)
                                jr = jr[ jr .∉ Ref(rmidx) ]  # keep those who are not to be deleted
                                if length(jr) > 0
                                    jr = minimum(jr)   # smallest of those
                                    # push!(insert_right, interp(m.c[id,iy,it].env[jr-1:jr], [ I.x ] ) )
                                    newright = MLine(m.c[id,iy,it].env.v[jr-1:jr])
                                    sortx!(newright)
                                    push!(insert_right, getv(interp(newright, [ I.x ] ))[1] )
                                # else
                                #     push!(insert_right,I)
                                end
                            end # all intersections

                            # remove illegal points from c
                            deleteat!(m.c[id,iy,it].env.v, rmidx)

                            # add new points in twice with a slight offset from left
                            # to preserve the ordering in x.
                            for ileft in 1:length(insert_left)

                                consx = getx(m.c[id,iy,it].env)
                                j = findfirst(consx .> insert_left[ileft].x)  # first point past new intersection
                                insert!(m.c[id,iy,it].env.v,j,insert_left[ileft])  # item is j-th index
                                insert!(m.c[id,iy,it].env.v,j+1,insert_right[ileft])
                            end


                        end
                    else   # if id==1
                        m.v[id,iy,it] = Envelope(vline)
                    end

                    # store the expected value at the lower boundary
                    # in a separate object
                    m.v[id,iy,it].vbound = ev[1]

                    # this creates the credit constrained region
                    prepend!(m.c[id,iy,it].env,[Point(m.avec[1],0.0)])
                    prepend!(m.v[id,iy,it].env,[Point(m.avec[1],ev[1])])
                    # if !issorted(m.c[id,iy,it].env)
                    #     println(isecs)
                    #     xx = getx(m.c[id,iy,it].env)
                    #     ii = findfirst( vcat(0,diff(xx)) .< 0 )
                    #     println(m.c[id,iy,it].env.v[ii-1:ii+1])
                        sortx!(m.c[id,iy,it].env)  # sort cons by default
                    # end
                    # @assert issorted(m.c[id,iy,it].env)
                    # sortx!(m.c[id,it].env)
                    # sortx!(m.v[id,it].env)
                    # prepend!(m.v[id,it].env,p.a_low,ev[1])
                    # do NOT prepend the value function with the special value from above.
                end # if last period
            end  # current id
        end  # iy
    end # it

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


function runf(;par=Dict())
    p = Param(par=par)
    p.beta = 1/p.R
    m = FModel(p)
    dc_EGM!(m,p)
    (m,p)
end

function rung(;par=Dict())
    p = Param(par=par)
    p.beta = 1/p.R
    m = GModel(p)
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

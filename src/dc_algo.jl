

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
                # initialize value function with bequest function
                m.v[id,it] = Envelope(MLine(m.avec, bequest(m.avec,p) ) )
            end
        else
            for id in 1:p.nD   # current period dchoice
                working = id==1  # working today is id=1
                # what's next period's cash on hand given you work/not TODAY?
                mm1 = m.m1[it][id]
                # clamp!(mm1, p.cfloor, Inf)  # dont actullay want to do that: want to allow for neg assets

                # next period consumption and values y-coords
                # for each d-choice
                cmat = fill(-Inf,p.nD,p.na*p.ny)
                vmat = fill(-Inf,p.nD,p.na*p.ny)

                # get future cons and values
                # for future nD compute consumption and value function
                for iid in 1:p.nD
                    c1 = interp(m.c[iid,it+1].env, mm1[:])
                    floory!(c1,p.cfloor)   # floor negative consumption
                    cmat[iid,:] = gety(c1)  # get y-values
                    vmat[iid,:] = vfun(iid,it+1,cmat[iid,:],mm1[:],m.v[iid,it+1],p)
                end

                # get ccp to be a worker
                pwork = working ? ccp(vmat,p) : p.delta*ones(size(vmat)[2])

                # get expected marginal utility of that consumption
                mu1 = reshape(pwork .* up(cmat[1,:],p) .+ (1.0 .- pwork) .* up(cmat[2,:],p),p.ny,p.na)

                # get expected marginal value of saving: RHS of euler equation
                # beta * R * E[ u'(c_{t+1}) ]
                # ywgt: weight on the income states
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
                    ev =  (1-p.delta)*m.ywgt' * reshape(vmat[2,:],p.ny,p.na)+p.delta*m.ywgt' * reshape(logsum(vmat,p),p.ny,p.na)
                end
                vline = MLine(m.avec .+ c0, u(c0,id==1,p) .+ p.beta * ev[:])

                # SECONDARY ENVELOPE COMPUTATION
                # ==============================

                if id==1   # only for workers
                    m.v[id,it], m.c[id,it] = do_secondary(vline,cline,working,ev[1],p)
                else   # if id==1
                    m.v[id,it] = Envelope(vline)
                    m.c[id,it] = Envelope(cline)
                end

                # store the expected value at the lower boundary
                # in a separate object
                m.v[id,it].vbound = ev[1]

                # this creates the credit constrained region
                prepend!(m.c[id,it].env,[Point(m.avec[1],0.0)])
                prepend!(m.v[id,it].env,[Point(m.avec[1],ev[1])])
                sortx!(m.c[id,it].env)  # sort cons by default
            end # current discrete choice
        end    # if final period
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
    # vtmp = fill(-Inf,p.nD,p.ny,p.na)

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
                    fill!(vmat,0.0)
                    fill!(ctmp,-Inf)
                    # fill!(vtmp,0.0)

                    for jy in 1:p.ny # future state: owner, renter, income, etc
                        pr = m.ywgt[iy,jy]  # transprob

                        for iid in 1:p.nD  # future dchoice
                            # only feasible choices at this state
                            # if renter, cannot sell etc

                            m1 = m.m1[it+1][iid][jy,:]  # state specific mvec
                            c1 = interp(m.c[iid,jy,it+1].env, m1) # C(d',y',m')
                            floory!(c1,p.cfloor)   # floor negative consumption
                            ctmp[iid,jy,:] = gety(c1)
                            # vtmp[iid,jy,:] = vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                            vmat[iid,:] += pr * vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                        end
                    end # end future state

                    # now get expectated value function conditional on iy: E[V(t+1,iid,y')|iy]
                    # vmat = integrate over second dimension
                    # vmat = dropdims( reduce(+, vtmp, dims = 2), dims = 2)

                    # get ccp of choices: P(d'|iy), pwork
                    pwork = working ? ccp(vmat,p) : zeros(size(vmat)[2])

                    # get y-expected MU of cons(t+1): uprime
                    up!(ctmp,p)

                    # integrate to get E[ u'(c(t+1,y')) | y]
                    # prepare
                    fill!(cmat,0.0)
                    for jd in 1:p.nD
                        for ja in 1:p.na
                            for jy in 1:p.ny
                                # cmat[jd,ja] += m.ywgt[iy,jy] * up(ctmp[jd,jy,ja],p)
                                cmat[jd,ja] += m.ywgt[iy,jy] * ctmp[jd,jy,ja]
                            end
                        end
                    end

                    # compute tomorrow's marginal utility
                    mu1 = pwork .* cmat[1,:] .+ (1.0 .- pwork) .* cmat[2,:] # 1,na

                    #RHS
                    RHS = p.beta * p.R * mu1

                    #optimal cons
                    c0 = iup(RHS,p)
                    # set optimal consumption function today. endo grid m and cons c0
                    cline = MLine(m.avec .+ c0, c0)
                    # store
                    # m.c[id,iy,it] = Envelope(cline)
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
                        m.v[id,iy,it], m.c[id,iy,it] = do_secondary(vline,cline,working,ev[1],p)
                    else   # if id==1
                        m.v[id,iy,it] = Envelope(vline)
                        m.c[id,iy,it] = Envelope(cline)
                    end

                    # store the expected value at the lower boundary
                    # in a separate object
                    m.v[id,iy,it].vbound = ev[1]


                    # this creates the credit constrained region
                    prepend!(m.c[id,iy,it].env,[Point(m.avec[1],0.0)])
                    prepend!(m.v[id,iy,it].env,[Point(m.avec[1],ev[1])])
                    sortx!(m.c[id,iy,it].env)  # sort cons by default
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

    if it == p.nT
        r[:] = bequest(m1,p)
    else
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
    end
    return r
end

function vfun(ufun::Function,it::Int,c1::Vector{Float64},m1::Vector{Float64},v::Envelope,p::Param)

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
        r[:] = ufun(c1) .+ p.beta * bound(v)
    elseif any(mask)
        r[mask] = ufun(c1[mask]) .+ p.beta * bound(v)
        # elsewhere
        r[.!mask] = gety(interp(v.env,m1[.!mask]))
    else
        r[:] = gety(interp(v.env,m1))
    end

    return r
end

"""
vfun used in simulation of state dependent model
"""
function vfun(ufun::Function,it::Int,c1::Vector{Float64},m1::Vector{Float64},vv::Vector{Envelope},p::Param)

    # each entry of c1,m1 and v represents an individual
    N = length(c1)
    @assert N == length(m1) == length(vv)
    # L = en.L[id]
    # v = en[id,iy,it]

    # computes v_{it}(m) = u(c) + beta v_{it+1}(m1)

    r = fill(NaN,N)
    critx = map(x -> getx(x.env)[2],vv)
    for (i,v) in enumerate(vv)
        if length(getx(v.env)) < 2
            error("need more than 2 points in envelope object")
        end
        if m1[i] < critx[i]  # credit constrained
        # if m1[i] < getx(v.env)[2]  # credit constrained
            r[i] = ufun(c1[i]) + p.beta * bound(v)
        else
            r[i] = gety(interp(v.env,[m1[i]]))[1]
        end
    end
    r
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
    # p.beta = 1/p.R
    m = FModel(p)
    dc_EGM!(m,p)
    (m,p)
end
function runfp(;par = Dict())
    m,p = runf(par = par)
    plot(m,p)
end



function rung(;par=Dict())
    p = Param(par=par)
    # p.beta = 1/p.R
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

function do_secondary(vline::MLine, cline::MLine, working::Bool, ev0::Float64, p::Param)
    minx = min_x(vline)
    if vline.v[1].x <= minx
        # normal case - no bending back behind first point
        # call secondary_envelope on vline
        vout = secondary_envelope(vline)
        cout = Envelope(cline)

    else
        # non-convex region lies inside credit constraint.
        # endogenous x grid bends back before the first x grid point.
        # println("minx = $minx")
        # println("first point = $(vline.v[1].x)")
        x0 = collect(range(minx,stop = vline.v[1].x,length = floor(Integer,p.na/10))) # some points to the left of first x point
        x0 = x0[1:end-1]
        c0 = copy(x0)
        if minx < 0
            c0[:] .= c0 .+ abs(minx) .+ p.cfloor
        end
        y0 = u(c0,working,p) .+ p.beta .* ev0   # use c0: positive consumption even with neg assets!

        prepend!(vline,convert(Point,x0,y0))
        prepend!(cline,convert(Point,x0,c0))  # cons policy in credit constrained is 45 degree line
        cout = Envelope(cline)
        vout = secondary_envelope(vline)
    end

    # now we have a cleaned value function in vout

    # if any points were removed from vline, now we need to clean the consumption function as well
    if length(getr(vout)) > 0

        # analyse policy function
        # =======================

        rmidx = getr(vout)  # indices removed from value function
        insert_left = Point[]
        insert_right = Point[]

        # remove illegal points from c
        # deleteat!(m.c[id,iy,it].env.v, rmidx)

        # insert new intersections into consumption function
        isecs = gets(vout)
        consx = getx(cout.env)

        # compute any new cons-points (intersections in v)
        # ---------------------------------------

        for isec in 1:length(isecs)
            I = isecs[isec]

            # interpolate from left
            jl = findall(consx .< I.x)
            jl = jl[ jl .∉ Ref(rmidx) ]  # keep those who are not to be deleted
            if length(jl) > 0
                jl = maximum(jl)  # biggest of those
                if jl < length(consx)
                    newleft = MLine(cout.env.v[jl:jl+1])
                    # @fediskhakov: who guarantees that jl+1 is not to be deleted?
                    # more generally: why do we not delete points before we
                    # insert the new intersections?

                    sortx!(newleft)
                    tmp = getv(interp(newleft, [ I.x ] ))[1]
                    # take this new point at a minimal left shift in x
                    push!(insert_left, Point(tmp.x - 1e3 * eps(), tmp.y) )
                # else
                #     println("trying to insert after last point")
                end
            # else
            #     push!(insert_left,I)
            end

            # interpolate from right
            jr = findall(consx .> I.x)
            jr = jr[ jr .∉ Ref(rmidx) ]  # keep those who are not to be deleted
            if length(jr) > 0
                jr = minimum(jr)   # smallest of those
                # push!(insert_right, interp(m.c[id,iy,it].env[jr-1:jr], [ I.x ] ) )
                newright = MLine(cout.env.v[jr-1:jr])
                sortx!(newright)
                push!(insert_right, getv(interp(newright, [ I.x ] ))[1] )
            # else
            #     push!(insert_right,I)
            end
        end # all intersections

        # remove illegal points from c
        deleteat!(cout.env.v, rmidx)

        # add new points in twice with a slight offset from left
        # to preserve the ordering in x.
        for ileft in 1:length(insert_left)

            consx = getx(cout.env)
            j = findfirst(consx .> insert_left[ileft].x)  # first point past new intersection
            insert!(cout.env.v,j,insert_left[ileft])  # item is j-th index
            insert!(cout.env.v,j+1,insert_right[ileft])
        end
    end
    (vout, cout)
end

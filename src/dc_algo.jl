


function minimal_EGM()
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))
    avec          = scaleGrid(p.a_low,p.a_high,p.na,logorder = 1)
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,xticks = 1:10)
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 1.0 .+ exp.(yvec).*p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
        c[it] = vcat(0, 1 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector. (p.na+1,1)
        m[it] = vcat(p.a_low, avec .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        plot!(pl,m[it],c[it],label="$it")
    end
    return (m,c,pl)
end



"""
    vfun(id::Int,iy::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Array{Envelope},p::Param)

Calculate the period `it`, `iy`-state, discrete choice `id`-specific value function. Avoids interpolation in credit constrained region by using the analytic form of the value function (no need to interpolate expected value function when on the lower bound of assets.)
"""
function vfun(id::Int,iy::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Array{Envelope,N} where N,p::Param)

    # L = en.L[id]
    v = en[id,iy,it]

    # computes v_{it}(m) = u(c) + beta v_{it+1}(m1)

    if length(getx(v.env)) < 2
        error("need more than 2 points in envelope object")
    end

    r = fill(NaN,size(m1))
    mask = m1.<getx(v.env)[2]
    mask = it==p.nT ? trues(size(mask)) : mask

    if all(mask)
        # in the credit constrained region:
        r[:] = u(c1,id==2,p) .+ p.beta * bound(v)
    elseif any(mask)
        r[mask] = u(c1[mask],id==2,p) .+ p.beta * bound(v)
        # elsewhere
        r[.!mask] = gety(interp(v.env,m1[.!mask]))
    else
        r[:] = gety(interp(v.env,m1))
    end

    return r
end

"""
    ccp(m::Model,p::Param)

Conditional Choice probability of working
"""
function ccp(x::Matrix,p::Param)
    #choice probability of the second row in 2-row matrix
    mx = maximum(x,dims = 1)
    mxx = x.-repeat(mx,size(x)[1],1)   # center values at max for numerical stability
    vec(exp.(mxx[2,:]./p.lambda)' ./ sum(exp.(mxx./p.lambda),dims = 1))
end

"""
    logsum(m::Model,p::Param)

Logsum of conditional values used in Expected value function.
"""
function logsum(x::Matrix,p::Param)
    mx = maximum(x,dims = 1)
    mxx = x .- repeat(mx,size(x)[1],1)
    mx .+ p.lambda * log.( sum(exp.(mxx./p.lambda), dims = 1) )
end


"""
    dcegm!(m::Model,p::Param)

Main body of the DC-EGM algorithm version
"""
function dc_EGM!(m::Model,p::Param)
    for it in p.nT:-1:1
        # println()
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
                    working = id==2  # working today is id=2
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
                            c1 = interp(m.c[iid,iy,it+1].env, mm1[:])
                            floory!(c1,p.cfloor)   # floor negative consumption
                            cmat[iid,:] = gety(c1)  # get y-values

                            vmat[iid,:] = vfun(iid,iy,it+1,cmat[iid,:],mm1[:],m.v,p)
                        end

                    else  # retirees have no discrete choice - absorbing state

                        iid = 1   # no work next period
                        mm1 = m.m1[it+1][iid]  # non-work wealth
                        c1 = interp(m.c[iid,iy,it+1].env, mm1[:])
                        floory!(c1,p.cfloor)   # floor negative consumption
                        cmat[iid,:] = gety(c1)

                        vmat[iid,:] = vfun(iid,iy,it+1,cmat[iid,:],mm1[:],m.v,p)

                    end

                    # get ccp to be a worker
                    pwork = working ? ccp(vmat,p) : zeros(size(vmat)[2])

                    # get marginal utility of that consumption
                    mu1 = reshape(pwork .* up(cmat[2,:],p) .+ (1.0 .- pwork) .* up(cmat[1,:],p),p.na,p.ny)
                    # println("mu1 = ")
                    # display(mu1[1:10,:])

                    # get expected marginal value of saving: RHS of euler equation
                    # beta * R * E[ u'(c_{t+1}) | iy ]
                    # need to integrate out Py here
                    RHS = p.beta * p.R * mu1 * m.ywgt[:,iy]
                    # println("RHS = $(RHS[1:10])")

                    # optimal consumption today: invert the RHS of euler equation
                    c0 = iup(RHS,p)

                    # set optimal consumption function today. endo grid m and cons c0
                    cline = MLine((m.avec[it]) .+ c0, c0)
                    # store
                    m.c[id,iy,it] = Envelope(cline)

                    # consumption function done.


                    # compute value function
                    # ----------------------
                    if working
                        ev = reshape(logsum(vmat,p),p.na,p.ny) * m.ywgt[:,iy]
                    else
                        ev = reshape(vmat[1,:],p.na,p.ny) * m.ywgt[:,iy]
                    end
                    vline = MLine((m.avec[it]) .+ c0, u(c0,id==2,p) .+ p.beta * ev)

                    # println(vline)

                    if any(isnan.(ev))
                        println("ev = ")
                        display(ev)
                    end



                    # vline and cline may have backward-bending regions: let's prune those
                    # SECONDARY ENVELOPE COMPUTATION

                    if id==2   # only for workers
                        minx = min_x(vline)
                        if minx < vline.v[1].x
                            # non-convex region lies inside credit constraint.
                            # endogenous x grid bends back before the first x grid point.
                            x0 = range(minx,stop = vline.v[1].x,length = floor(p.na/10)) # some points to the left of first x point
                            x0 = x0[1:end-1]
                            y0 = u(x0,working,p) + p.beta * ev[1]
                            prepend!(vline,[Point(x0,y0)])
                            prepend!(cline,[Point(x0,y0)])  # cons policy in credit constrained is 45 degree line
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


function run2()
    p = Param()
    m = Model(p)
    dc_EGM!(m,p)
    (m,p)
    # plot(m,id=2,xlim=(0,4),ylim=(-100,-20))
    # plot!(m,id=1,linestyle=:dash,label="")
end
function runit()
    p = Param()
    m = Model(p)
    dc_EGM!(m,p)
    return (m,p)
end
function pp()
    m,p = runit()
    plot(m.v[1,1].env)
end

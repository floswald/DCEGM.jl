
# status: WIP
# still does not allow to come back from bankruptcy punishment
# this means that final period bankrtuptcy punishment is very important
function bk!(m::BModel,p::Param)

    cmat = fill(-Inf,p.nD,p.na)
    cmat2 = fill(-Inf,p.nD,p.na)
    vmat = fill(-Inf,p.nD,p.na)
    ctmp = fill(-Inf,p.nD,p.ny,p.na)
    # vtmp = fill(-Inf,p.nD,p.ny,p.na)

    for it in p.nT:-1:1
        # println(it)
        for iy in 1:p.ny  # current state
            for id in 1:p.nD  # current dchoice
                # println("id = $id")
                filer = id==2  # filing is id=2

                if it==p.nT
                    # final period: consume everyting.
                    if id == 1
                        m.c[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 , p.a_high)) )
                        m.v[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 , NaN)) )
                    else
                        m.c[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 , p.a_high)) )
                        m.v[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 + p.alphaT  , NaN)) )
                        m.v[id,iy,it].vbound = p.alphaT
                    end
                    # initialize value function with vf(1) = 0
                    # m.v[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(id == 1 ? 0.0 : p.alphaT,p.alphaT)) )

                    # with bk flag on:
                    m.cbk[iy,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(0.0,p.a_high)) )
                    # initiz,iae value function with vf(1) = 0
                    m.vbk[iy,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(p.alphaT,NaN)) )
                    m.vbk[iy,it].vbound = p.alphaT
                else
                    # next period consumption and values y-coords
                    # for each d-choice
                    # reset all value matrices
                    fill!(vmat,0.0)
                    fill!(ctmp,-Inf)
                    fill!(cmat,0.0)
                    fill!(cmat2,0.0)

                    # fill!(vtmp,0.0)#


                    # you are the bank
                    # you know how pairs (ia,iy) map into bankruptcies in any given period
                    # like you observe people in t+1 and you can now know that a certain asset state ia together with iy
                    # corresponds to "bankruptcy".
                    # in period t, you can compute an expectation over this.
                    prnobk_at_ia = zeros(p.na)
                    for jy in 1:p.ny # future state: owner, renter, income, etc
                        # E[ m.v[1,jy,it+1] vs m.v[2,jy,it+1] | iy ]
                        # m.v is defined on cash on hand M = a + y, so I can evaluate it at any M
                        # say M = iy + ia
                        pr = m.ywgt[iy,jy]  # transprob
                        # println("v1 = $(gety(interp(m.v[1,jy,it+1].env, income(it+1,p,m.yvec[jy]) .+ m.avec )))")
                        if it+1 == p.nT
                            I = trues(p.na)
                        else
                            may = income(it+1,p,m.yvec[jy]) .+ m.avec
                            b = [gety(interp(m.v[1,jy,it+1].env, may )) gety(interp(m.v[2,jy,it+1].env, may))]
                            b[isnan.(b)] .= -Inf
                            I = b[:,1] .> b[:,2]
                        end
                        prnobk_at_ia[:] += pr * I
                    end

                    # risk adjusted return times asset grid
                    m.rgrid[id,iy,it] = clamp!([(p.R / prnobk_at_ia[ia]) for ia in 1:p.na], p.R, p.maxinterest)
                    # rgrid = p.R
                    # println("rgrid = $rgrid")

                    ragrid = m.rgrid[id,iy,it] .* m.avec

                    if !filer
                        fill!(cmat,0.0)

                        for jy in 1:p.ny # future state: owner, renter, income, etc
                            pr = m.ywgt[iy,jy]  # transprob

                            for iid in 1:p.nD  # future dchoice
                                # only feasible choices at this state
                                # if renter, cannot sell etc
                                if iid == 1
                                    # cash on hand : R (a[ia])*a[ia] + y(t+1)
                                    # R[ia] = (1+r) / pr_no_bk_iy[ia]
                                    # cash on hand:
                                    m1 = ragrid .+ income(it+1,p,m.yvec[jy])
                                    # m1 = m.m1[it+1][iid][jy,:]  # state specific mvec
                                    c1 = interp(m.c[iid,jy,it+1].env, m1) # C(d',y',m')
                                    floory!(c1,p.cfloor)   # floor negative consumption
                                    ctmp[iid,jy,:] = gety(c1)
                                    # vtmp[iid,jy,:] = vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                                    vmat[iid,:] += pr * vfun(x->u(x,filer,p),it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                                    cmat[iid,:] += pr * up(ctmp[iid,jy,:],p)
                                else
                                    # you are a filer: no savings and analytic value function
                                    ctmp[iid,jy,:] .= income(it+1,p,repeat([m.yvec[jy]],p.na))
                                    vmat[iid,:] += pr * (u(ctmp[iid,jy,:],true,p) .+ p.beta * bound(m.vbk[jy,it+1]))
                                    cmat[iid,:] += pr * up(ctmp[iid,jy,:],p)
                                    vmat[iid,m.iazero:end] .= -Inf   # filing only with negative assets!
                                    cmat[iid,m.iazero:end] .= 0.0   # filing only with negative assets!
                                end
                            end
                        end # end future state

                        # Discrete choice setup next period
                        # =================================

                        # get ccp of filing next period at m1 : P(d'|iy), pfile
                        pnofile = ccp(vmat,p)   # prob of filing TOMORROW from today's perspective
                        # m.pnofile[iy,it,:] = pnofile   # prob of filing at state ia next period, given iy in period it.

                        # compute tomorrow's marginal utility
                        mu1 = pnofile .* cmat[1,:] .+ (1.0 .- pnofile) .* cmat[2,:] # 1,na

                        #RHS
                        # TODO R here will be computed above as the expected return of the bank
                        # in each discrete choice
                        RHS = p.beta * p.R * mu1

                        #optimal cons
                        c0 = iup(RHS,p)
                        # set optimal consumption function today. endo grid m and cons c0
                        cline = MLine(m.avec .+ c0, c0)
                        # consumption function done.

                        # compute value function
                        # ----------------------
                        ev =  logsum(vmat,p)  # E[ max(V(1,y'),V(2,y')) | y]
                        vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev[:])
                        minv = gety(vline)[1]

                        # store the expected value at the lower boundary
                        # in a separate object
                        minv = ev[1] < minv ? ev[1] : minv

                        # SECONDARY ENVELOPE COMPUTATION
                        # ==============================
                        m.v[id,iy,it], m.c[id,iy,it] = do_secondary(vline,cline,filer,minv,p)
                        m.v[id,iy,it].vbound = minv


                        # this creates the credit constrained region
                        prepend!(m.c[id,iy,it].env,[Point(m.avec[1],0.0)])
                        prepend!(m.v[id,iy,it].env,[Point(m.avec[1],m.v[id,iy,it].vbound)])
                        sortx!(m.c[id,iy,it].env)  # sort cons by default

                    else  # if you are filing today (or you have filed previously)

                        iid = 2 # hard code to 2 - we only will use the second row of vmat here
                        # there is no discrete choice for you next period.

                        # we compute the value and savings function of being in bk state first.
                        for jy in 1:p.ny # future state:
                            pr = m.ywgt[iy,jy]  # transprob
                            m1 = m.m1[it+1][jy,:]  # cash on hand in bk state
                            c1 = interp(m.cbk[jy,it+1].env, m1) # C(d',y',m')
                            floory!(c1,p.cfloor)   # floor negative consumption
                            ctmp[iid,jy,:] = gety(c1)
                            # vtmp[iid,jy,:] = vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                            vmat[iid,:] += (pr * vfun(x->u(x,filer,p),it+1,ctmp[iid,jy,:],m1,m.vbk[jy,it+1],p))
                            cmat[iid,:] += (pr * up(ctmp[iid,jy,:],p))  # E[ c(id',y') | y]

                            # # now get the same for non-bk state
                            # c1 = interp(m.c[1,jy,it+1].env, m1) # C(d',y',m')
                            # floory!(c1,p.cfloor)   # floor negative consumption
                            # ctmp[iid,jy,:] = gety(c1)
                            # # vtmp[iid,jy,:] = vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                            # vmat[iid,:] += (pr * vfun(x->u(x,filer,p),it+1,ctmp[iid,jy,:],m1,m.vbk[jy,it+1],p))
                            # cmat[iid,:] += (pr * up(ctmp[iid,jy,:],p))  # E[ c(id',y') | y]

                        end  # future states

                        mu1 = cmat[iid,:]
                        #RHS
                        RHS = p.beta * p.R * mu1
                        #optimal cons
                        c0 = iup(RHS,p)
                        # set optimal consumption function today. endo grid m and cons c0
                        cline = MLine(m.aposvec .+ c0, c0)
                        # consumption function done.
                        ev = vmat[iid,:]
                        vline = MLine(m.aposvec .+ c0, u(c0,id==2,p) .+ p.beta * ev[:])

                        m.cbk[iy,it] = Envelope(cline)
                        m.vbk[iy,it] = Envelope(vline)
                        m.vbk[iy,it].vbound = ev[1]

                        # this creates the credit constrained region
                        prepend!(m.cbk[iy,it].env,[Point(m.aposvec[1],0.0)])
                        prepend!(m.vbk[iy,it].env,[Point(m.aposvec[1],ev[1])])
                        sortx!(m.cbk[iy,it].env)  # sort cons by default
                        # end value in bk state

                        # value and cons function if filing today
                        c0 = income(it,p,repeat([m.yvec[iy]],p.na))
                        c01 = copy(c0)
                        c01[m.iazero:end] .= NaN
                        cline = MLine(m.avec .+ c0, c01)
                        vline = MLine(m.avec .+ c0, u(c01, id == 2, p) .+ p.beta * ev[1])
                        m.c[id,iy,it] = Envelope(cline)
                        m.v[id,iy,it] = Envelope(vline)
                        prepend!(m.c[id,iy,it].env,[Point(m.avec[1],c01[1])])
                        prepend!(m.v[id,iy,it].env,[Point(m.avec[1],gety(vline)[1])])

                        m.v[id,iy,it].vbound = NaN   # does not exist - you don't save when filing.
                    end
                end # if last period
            end  # current id
        end  # iy
    end # it
end


function runbk(; par = Dict(:nT => 50, :a_low => -5.0,:a_lowT => -5.0,:na =>101, :alpha => 0.0,  :alphaT => 0.0, :lambda => 0.5))
    p = Param(par = par)
    m = BModel(p)
    bk!(m,p)
    m,p
end

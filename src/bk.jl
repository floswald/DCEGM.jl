

function bk!(m::BModel,p::Param)

    cmat = fill(-Inf,p.nD,p.na)
    vmat = fill(-Inf,p.nD,p.na)
    ctmp = fill(-Inf,p.nD,p.ny,p.na)
    # vtmp = fill(-Inf,p.nD,p.ny,p.na)

    for it in p.nT:-1:1
        println(it)
        for iy in 1:p.ny  # current state
            for id in 1:p.nD  # current dchoice
                println("id = $id")
                filer = id==2  # filing is id=2

                if it==p.nT
                    # final period: consume everyting.
                    m.c[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(id == 1 ? 0.0 : NaN,p.a_high)) )
                    # initialize value function with vf(1) = 0
                    m.v[id,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(id == 1 ? 0.0 : NaN,NaN)) )
                    # note that 0.0 as first value of the vfun is not innocuous here!

                    # with bk flag on:
                    m.cbk[iy,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(0.0,p.a_high)) )
                    # initiz,iae value function with vf(1) = 0
                    m.vbk[iy,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(-p.alpha,NaN)) )
                    m.vbk[iy,it].vbound = -p.alpha
                else
                    # next period consumption and values y-coords
                    # for each d-choice
                    # reset all value matrices
                    fill!(vmat,0.0)
                    fill!(ctmp,-Inf)
                    fill!(cmat,0.0)

                    # fill!(vtmp,0.0)

                    if !filer
                        for jy in 1:p.ny # future state: owner, renter, income, etc
                            pr = m.ywgt[iy,jy]  # transprob

                            for iid in 1:p.nD  # future dchoice
                                # only feasible choices at this state
                                # if renter, cannot sell etc
                                if iid == 1
                                    m1 = m.m1[it+1][iid][jy,:]  # state specific mvec
                                    c1 = interp(m.c[iid,jy,it+1].env, m1) # C(d',y',m')
                                    floory!(c1,p.cfloor)   # floor negative consumption
                                    ctmp[iid,jy,:] = gety(c1)
                                    # vtmp[iid,jy,:] = vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                                    vmat[iid,:] += pr * vfun(x->u(x,filer,p),it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                                else
                                    # you are a filer: no savings and analytic value function
                                    ctmp[iid,jy,:] .= income(it+1,p,repeat([m.yvec[jy]],p.na))
                                    vmat[iid,:] += pr * (u(ctmp[iid,jy,:],true,p) .+ p.beta * bound(m.vbk[jy,it+1]))
                                    vmat[iid,m.iazero:end] .= -Inf   # filing only with negative assets!
                                    ctmp[iid,jy,m.iazero:end] .= -Inf   # filing only with negative assets!
                                end
                            end
                        end # end future state

                        # Discrete choice setup next period
                        # =================================

                        # get ccp of filing next period: P(d'|iy), pfile
                        pnofile = ccp(vmat,p)

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
                        mu1 = pnofile .* cmat[1,:] .+ (1.0 .- pnofile) .* cmat[2,:] # 1,na

                        #RHS
                        RHS = p.beta * p.R * mu1

                        #optimal cons
                        c0 = iup(RHS,p)
                        # set optimal consumption function today. endo grid m and cons c0
                        cline = MLine(m.avec .+ c0, c0)
                        # consumption function done.

                        # compute value function
                        # ----------------------
                        ev =  logsum(vmat,p)
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
                            m1 = m.m1[it+1][iid][jy,:]  # cash on hand in bk state
                            c1 = interp(m.cbk[jy,it+1].env, m1) # C(d',y',m')
                            floory!(c1,p.cfloor)   # floor negative consumption
                            ctmp[iid,jy,:] = gety(c1)
                            # vtmp[iid,jy,:] = vfun(iid,it+1,ctmp[iid,jy,:],m1,m.v[iid,jy,it+1],p)
                            vmat[iid,:] += pr * vfun(x->u(x,filer,p),it+1,ctmp[iid,jy,:],m1,m.vbk[jy,it+1],p)
                            cmat[iid,:] += pr * up(ctmp[iid,jy,:],p)
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


function runbk()
    p = Param(par = Dict(:a_low => -5.0,:na =>101, :alpha => 0.1))
    m = BModel(p)
    bk!(m,p)
    m,p
end



"""
    bkdchoice!(ctmp::Array,vmat::Array,p::Param,m::Model,it::Int,jy::Int,pr::Float64)

fills preallocated ctmp and vmat arrays with state-specific future consumption and value functions _conditional_
on future discrete choices taken.

* `ctmp`: next period consumption in state `jy` and dchoice `iid`
* `vmat`: next period value in dchoice `iid`. Expectation wrt `jy` integrated out.
"""
function bkdchoice!(ctmp::Array,vmat::Array,p::Param,m::Model,m1::Vector{Float64},it::Int,jy::Int,pr::Float64)
    for iid in 1:p.nD  # future dchoice
        # only feasible choices at this state
        # if renter, cannot sell etc
        # TODO improve hard coding of iflag here
        iflag = 1
        if iid == 1
            c1 = interp(m.c[iid,iflag,jy,it+1].env, m1) # C(d',y',m')  iflag = 1
            floory!(c1,p.cfloor)   # floor negative consumption
            ctmp[iid,jy,:] = gety(c1)
            # filer = false
            vmat[iid,:] += pr * vfun(x->u(x,false,p),it+1,ctmp[iid,jy,:],m1,m.v[iid,iflag,jy,it+1],p)
        else
            # you are filing next period : no savings and analytic value function
            ctmp[iid,jy,:] .= income(it+1,p,repeat([m.yvec[jy]],p.na))
            vmat[iid,:] += pr * (u(ctmp[iid,jy,:],true,p) .+ p.beta * bound(m.v[1,2,jy,it+1]))
            vmat[iid,m.iazero:end] .= -Inf   # filing only with negative assets!
            ctmp[iid,jy,m.iazero:end] .= -Inf   # filing only with negative assets!
        end
    end
end


function bk!(m::BModel,p::Param)
    # TODO need to make distinction between `flag` and `dchoice`
    # some states (flag) allow some dchoices while others dont

    cmat = fill(-Inf,p.nD,p.na)
    vmat = fill(-Inf,p.nD,p.na)
    ctmp = fill(-Inf,p.nD,p.ny,p.na)
    # vtmp = fill(-Inf,p.nD,p.ny,p.na)

    cvec = fill(-Inf,p.na)
    vvec = fill(-Inf,p.na)
    vvec1 = fill(-Inf,p.na) # future value of coming back from punishment
    ctmp1 = fill(-Inf,p.ny,p.na)

    for it in p.nT:-1:1
        println("it = $it")
        for iy in 1:p.ny  # current state
            println("iy = $iy")
            for iflag in 1:2  # iflag = 1: not in bankruptcy state, iflag = 2: in bankruptcy state
                println("iflag = $iflag")
                if it==p.nT
                    # final period: consume everyting.
                    # if no flag, can choose to file (but get bad value)
                    if iflag == 1
                        for id in 1:p.nD
                            if id == 1
                                m.c[id,iflag,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 ,p.a_high)) )
                                m.v[id,iflag,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 ,NaN)) )
                            else
                                # if you file in last period, you get alphaT less utilty at each cash-on-hand value.
                                # that prevents anyone from choosing filing in last period.
                                m.c[id,iflag,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 ,p.a_high)) )
                                m.v[id,iflag,iy,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat( 0.0 + p.alphaT  ,NaN)) )
                                m.v[id,iflag,iy,it].vbound = p.alphaT
                            end
                        end
                    else
                        # if flag, only one option (id = 1), and get worse value.
                        m.c[1,iflag,iy,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(0.0,p.a_high)) )
                        # initiz,iae value function with vf(1) = 0
                        m.v[1,iflag,iy,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(-p.alpha,NaN)) )
                        m.v[1,iflag,iy,it].vbound = p.alphaT
                    end
                else  # preceding periods



                    if iflag == 1  # no flag: can choose to file or not to file today.

                        # for each d-choice
                        # reset all value matrices
                        fill!(vmat,0.0)
                        fill!(ctmp,-Inf)
                        fill!(cmat,0.0)

                        for id in 1:p.nD

                            if id == 1  # not filing today
                                # compute next period values and consumption functions
                                for jy in 1:p.ny # future state: owner, renter, income, etc
                                    pr = m.ywgt[iy,jy]  # transprob
                                    bkdchoice!(ctmp,vmat,p,m,m.m1[it+1][iflag][jy,:],it,jy,pr)   # TODO this should have iid instead of iflag!
                                end # end future state

                                # Probability of making Discrete choices next period
                                # ==================================================

                                # get ccp of filing next period: P(d'|iy), pfile
                                if any(isnan.(vmat))
                                    vmat[isnan.(vmat)] .= -Inf
                                end
                                pnofile = ccp(vmat,p)



                                # get y-expected MU of cons(t+1): uprime
                                up!(ctmp,p)
                                if any(isnan.(ctmp))
                                    println(ctmp)
                                    error()
                                end

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

                                ev =  logsum(vmat,p)  # expected value function (integration over EV1 shocks)
                                vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev[:])
                                minv = gety(vline)[1]

                                # store the expected value at the lower boundary
                                # in a separate object
                                minv = ev[1] < minv ? ev[1] : minv
                                if any(isnan.(ev))
                                    println(ev)
                                    error()
                                end

                                # SECONDARY ENVELOPE COMPUTATION
                                # ==============================
                                m.v[id,iflag,iy,it], m.c[id,iflag,iy,it] = do_secondary(vline,cline,false,minv,p)  # `false` not filing
                                m.v[id,iflag,iy,it].vbound = minv
                                m.ccp[id,iflag,iy,it] = Envelope(MLine(m.avec .+ c0, pnofile))
                                m.ev[id,iflag,iy,it] = Envelope(MLine(m.avec .+ c0, ev[:]))

                                # this creates the credit constrained region
                                prepend!(m.c[id,iflag,iy,it].env,[Point(m.avec[1],0.0)])
                                prepend!(m.v[id,iflag,iy,it].env,[Point(m.avec[1],m.v[id,iflag,iy,it].vbound)])
                                sortx!(m.c[id,iflag,iy,it].env)  # sort cons by default

                            # else  # filing for bk today
                            #     # no savings choice, consume income.
                            #     c0 = income(it,p,repeat([m.yvec[iy]],p.na))
                            #     c01 = copy(c0)
                            #     c01[m.iazero:end] .= NaN
                            #     cline = MLine(m.avec .+ c0, c01)
                            #     vline = MLine(m.avec .+ c0, u(c01, true, p) .+ p.beta * ev[1])  # TODO wrong ev!! needs m.v[id,iflag,iy,it].vbound
                            #     m.c[id,iflag,iy,it] = Envelope(cline)
                            #     m.v[id,iflag,iy,it] = Envelope(vline)
                            #     prepend!(m.c[id,iflag,iy,it].env,[Point(m.avec[1],c01[1])])
                            #     prepend!(m.v[id,iflag,iy,it].env,[Point(m.avec[1],gety(vline)[1])])
                            #     m.v[id,iflag,iy,it].vbound = NaN   # does not exist - you don't save when filing.

                            end

                        end # for all id when iflag == 1

                    else  # iflag == 2
                        # no discrete choices here
                        # so just will use id = 1 all the time
                        id = 1
                        iid = 1

                        # reset all value vectors
                        fill!(vvec,0.0)
                        fill!(vvec1,0.0)
                        fill!(ctmp1,-Inf)
                        fill!(cmat,0.0)

                        # we compute the next period value and savings function of remaining in bk state.
                        for jy in 1:p.ny # future y-state:
                            pr = m.ywgt[iy,jy]  # transprob
                            m1 = m.m1[it+1][iflag][jy,:]  # cash on hand in bk state
                            c1 = interp(m.c[id,iflag,jy,it+1].env, m1) # Consumption function in bk state (iflag = 2)
                            floory!(c1,p.cfloor)   # floor negative consumption
                            ctmp1[jy,:] = gety(c1)

                            # building future value functions
                            # `true` you are in punishment
                            # future value fun to interpolate is id=1, iflag=2 (i.e. the bk-state one)
                            vvec[:] += (pr * vfun(x->u(x,true,p),it+1,ctmp1[jy,:],m1,m.v[id,iflag,jy,it+1],p))

                            # but we also need the value of coming back to non bk-state, i.e. iflag = 1!
                            # interpolated at next period's cash on hand, which is just positive assets + income
                            vvec1[:] += (pr * gety(interp(m.v[1,1,iy,it+1].env,m1)))
                            cvec[:] += (pr * up(ctmp1[jy,:],p))
                        end  # future states

                        # TODO store vmat from above as m.ev[id,iflag,iy,it+1] = E[ V(d,it+1) | iy]
                        # and ccp. then can compute ev = ccp * vmat easier


                        #RHS
                        RHS = p.beta * p.R * cvec
                        #optimal cons
                        c0 = iup(RHS,p)
                        # set optimal consumption function today. endo grid m and cons c0
                        cline = MLine(m.aposvec .+ c0, c0)
                        # consumption function done.

                        # expected value function: there is random exit from the state with p.delta
                        ev = (1 - p.delta) * vvec + p.delta * vvec1
                        vline = MLine(m.aposvec .+ c0, u(c0,true,p) .+ p.beta * ev[:])

                        m.c[id,iflag,iy,it] = Envelope(cline)
                        m.v[id,iflag,iy,it] = Envelope(vline)
                        m.v[id,iflag,iy,it].vbound = ev[1]

                        # this creates the credit constrained region
                        prepend!(m.c[id,iflag,iy,it].env,[Point(m.aposvec[1],0.0)])
                        prepend!(m.v[id,iflag,iy,it].env,[Point(m.aposvec[1],ev[1])])
                        sortx!(m.c[id,iflag,iy,it].env)  # sort cons by default


                        # finish flag = 1, id =2 choice now
                        # no savings choice, consume income.
                        c0 = income(it,p,repeat([m.yvec[iy]],p.na))
                        c01 = copy(c0)
                        c01[m.iazero:end] .= NaN
                        cline = MLine(m.avec .+ c0, c01)
                        vline = MLine(m.avec .+ c0, u(c01, true, p) .+ p.beta * m.v[id,iflag,iy,it].vbound)
                        m.c[2,1,iy,it] = Envelope(cline)
                        m.v[2,1,iy,it] = Envelope(vline)
                        prepend!(m.c[2,1,iy,it].env,[Point(m.avec[1],c01[1])])
                        prepend!(m.v[2,1,iy,it].env,[Point(m.avec[1],gety(vline)[1])])
                        m.v[2,1,iy,it].vbound = NaN   # does not exist - you don't save when filing.


                    end  # end iflag == 2
                end # if last period
            end  # end iflag
        end  # iy
    end # it
end


function runbk(; par = Dict(:nT => 50, :a_low => -5.0,:a_lowT => 0.0,:na =>101, :alpha => 0.0,  :alphaT => 0.0, :lambda => 0.5))
    p = Param(par = par)
    m = BModel(p)
    bk!(m,p)
    m,p
end

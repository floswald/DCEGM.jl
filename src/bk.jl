

function bk!(m::BModel,p::Param)
    cmat = fill(-Inf,p.nD,p.na*p.ny)
    vmat = fill(-Inf,p.nD,p.na*p.ny)
    vtmp = fill(-Inf,p.nD,p.ny,p.na)
    for it in p.nT:-1:1
        # println(it)
        # @info("period = $it")

        if it==p.nT
            m.c[1,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,p.a_high)) )
            # initialize value function with vf(1) = 0
            m.v[1,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,NaN)) )
            m.c[2,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(0.0,p.a_high)) )  # this cons function will never be used.
            m.v[2,it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(NaN,NaN)) )   # cannot file in lastperiod

            # with bk flag on:
            m.cbk[it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(0.0,p.a_high)) )
            # initialize value function with vf(1) = 0
            m.vbk[it] = Envelope(MLine(vcat(0.0,p.a_high),vcat(0.0,NaN)) )

        else
            # bkflag off

                # choose file or not today
                id = 1
                    # if not file today, can choose file or not tomorrow
                    # if id == 1
                        # get future cons and values
                        iid = 1   # no file tomorrow
                        mm1 = m.m1[it][id]  # next period's cash on hand
                        c1 = interp(m.c[iid,it+1].env, mm1[:])
                        floory!(c1,p.cfloor)   # floor negative consumption
                        cmat[iid,:] = gety(c1)  # get y-values
                        vmat[iid,:] = vfun( x->u(x,false,p), it+1,cmat[iid,:],mm1[:],m.v[iid,it+1],p)  # false: not bankruptcy


                        iid = 2   # file tomorrow
                        c1 = interp(m.c[iid,it+1].env, mm1[:])    # filer's consumption function
                        floory!(c1,p.cfloor)   # floor negative consumption
                        cmat[iid,:] = gety(c1)  # get y-values
                        vtmp = reshape(vfun(x->u(x,true,p),it+1,cmat[iid,:],mm1[:],m.vbk[it+1],p),p.ny,p.na)  # next period value function for filer
                        vtmp[:,m.iazero:end] .= -Inf  # invalid choice for positive assets
                        vmat[iid,:] = vtmp[:]  # stick back into matrix

                        # vmat[iid,:] = u(cmat[iid,:],iid==2,p) .+ p.beta * vbound(m.vbk[it+1]))   # note the different future value here!
                        # vtmp = u(reshape(cmat[iid,:],p.ny,p.na),iid==2,p) .+ p.beta * vbound(m.vbk[it+1]))
                        # vtmp = u(reshape(cmat[iid,:],p.ny,p.na),iid==2,p) .+ p.beta * vbound(m.vbk[it+1]))
                        # now set values corresponding to a >= 0 to -Inf

                        # get ccp of those choices
                        pnobk = ccp(vmat,p)

                        # get expected marginal utility of that consumption
                        mu1 = reshape(pnobk .* up(cmat[1,:],p) .+ (1.0 .- pnobk) .* up(cmat[2,:],p),p.ny,p.na)

                        # get expected marginal value of saving: RHS of euler equation
                        # beta * R * E[ u'(c_{t+1}) ]
                        RHS = p.beta * p.R *  m.ywgt' * mu1

                        # optimal consumption today: invert the RHS of euler equation
                        c0 = iup(Array(RHS)[:],p)

                        # set optimal consumption function today. endo grid m and cons c0
                        cline = MLine(m.avec .+ c0, c0)
                        # store
                        m.c[id,it] = Envelope(cline)

                        # value function
                        ev =  m.ywgt' * reshape(logsum(vmat,p),p.ny,p.na)
                        vline = MLine(m.avec .+ c0, u(c0,true,p) .+ p.beta * ev[:])

                        # m.v[id,it], m.c[id,it] = do_secondary(vline,cline,id == 2,ev[1],p)  # id == 2 -> false (not in bk state)
                        m.v[id,it], m.c[id,it] = (Envelope(vline),Envelope(cline)) # id == 2 -> false (not in bk state)

                        m.v[id,it].vbound = ev[1]
                        prepend!(m.c[id,it].env,[Point(m.avec[1],0.0)])
                        prepend!(m.v[id,it].env,[Point(m.avec[1],ev[1])])
                        sortx!(m.c[id,it].env)  # sort cons by default



                id = 2
                        # else if file today, no dchoice tomorrow
                        mm1 = m.m1[it][id]  # next period's cash on hand: no negative assets possible
                        c1 = interp(m.cbk[it+1].env, mm1[:]) # correct cons functino
                        floory!(c1,p.cfloor)   # floor negative consumption
                        mu1 = reshape(up(gety(c1),p),p.ny,p.na)
                        RHS = p.beta * p.R *  m.ywgt' * mu1
                        # optimal consumption today: invert the RHS of euler equation
                        c0 = iup(Array(RHS)[:],p)
                        # set optimal consumption function today. endo grid m and cons c0
                        cline = MLine(m.avec .+ c0, c0)
                        # store
                        m.c[id,it] = Envelope(cline)

                        # next period value function if you file today : vtmp


                        # ev_file =  m.ywgt' * reshape(vtmp,p.ny,p.na)  #
                        ev_file =  m.ywgt' * reshape(vmat[2,:],p.ny,p.na)  #
                        vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev_file[:])   # is -Inf for a>0
                        m.v[id,it] = Envelope(vline)
                        m.v[id,it].vbound = ev_file[1]

                        prepend!(m.c[id,it].env,[Point(m.avec[1],0.0)])
                        prepend!(m.v[id,it].env,[Point(m.avec[1],ev_file[1])])
                        sortx!(m.c[id,it].env)  # sort cons by default


                        # at the same time store as cons function when bk flag is on
                        cline = MLine(m.aposvec .+ c0, c0)
                        # store
                        m.cbk[it] = Envelope(cline)

                        vfunbk = vfun( x->u(x,true,p),it+1,gety(c1),mm1[:],m.vbk[it+1],p)   # always pay the utility cost while in bk state
                        evbk =  m.ywgt' * reshape(vfunbk,p.ny,p.na)
                        vlinebk = MLine(m.aposvec .+ c0, u(c0,false,p) .+ p.beta * evbk[:])   # is -Inf for a>0
                        m.vbk[it] = Envelope(vlinebk)
                        m.vbk[it].vbound = evbk[1]

                        prepend!(m.cbk[it].env,[Point(m.aposvec[1],0.0)])
                        prepend!(m.vbk[it].env,[Point(m.aposvec[1],evbk[1])])
                        sortx!(m.cbk[it].env)  # sort cons by default



                  # end today's dchoice if bkflag is off

                # if you file for bk today.


                # compute value functions
                # ----------------------
                # if notbk
                #     ev =  m.ywgt' * reshape(logsum(vmat,p),p.ny,p.na)
                #     vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev[:])
                # else
                #     ev =  m.ywgt' * reshape(vmat[2,:],p.ny,p.na)
                #     vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev[:])
                # end

                # SECONDARY ENVELOPE COMPUTATION
                # ==============================

                # if id==1   # only if not bk
                #     m.v[id,it], m.c[id,it] = do_secondary(vline,cline,working,ev[1],p)
                # else   # if bk
                #     m.v[id,it] = Envelope(vline)
                #     m.c[id,it] = Envelope(cline)
                # end

                # store the expected value at the lower boundary
                # in a separate object
                # m.v[id,it].vbound = ev[1]

                # this creates the credit constrained region
                # prepend!(m.c[id,it].env,[Point(m.avec[1],0.0)])
                # prepend!(m.v[id,it].env,[Point(m.avec[1],ev[1])])
                # sortx!(m.c[id,it].env)  # sort cons by default
            # end # current discrete choice
        end    # if final period
    end     # loop over time
end

function runbk()
    p = Param(par = Dict(:a_low => -5.0))
    m = BModel(p)
    bk!(m,p)
    m,p
end

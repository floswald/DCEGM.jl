


function minimal_EGM(;dplot=false)
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    if dplot
        plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false)
    end
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 1.0 .+ exp.(yvec).*p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)  
        c1[c1.<0] .= p.cfloor
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
        c[it] = vcat(0, 1 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector. (p.na+1,1)
        m[it] = vcat(p.a_low, avec .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        if dplot
            plot!(m[it],c[it],label="$it")
        end
    end
    if dplot
        gui()
    end

    return (m,c)
end



"""
    eeee(id::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Matrix{Envelope},p::Param)

Calculate the period `it`, discrete choice `id`-specific value function. Avoids interpolation in credit constrained region by using the analytic form of the value function (no need to interpolate expected value function when on the lower bound of assets.)
"""
function vfun(id::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Matrix{Envelope},p::Param)

    # L = en.L[id]
    v = en[id,it]

    # computes v_{it}(m) = u(c) + beta v_{it+1}(m1)

    if length(getx(v)) < 2
        error("need more than 2 points in envelope object")
    end

    r = fill(NaN,size(m1))
    mask = m1.<getx(v)[2]
    mask = it==p.nT ? trues(size(mask)) : mask

    if all(mask)
        # in the credit constrained region:
        r[:] = u(c1,id==2,p) .+ p.beta * bound(v)
    elseif any(mask)
        r[mask] = u(c1[mask],id==2,p) .+ p.beta * bound(v)
        # elsewhere
        r[.!mask] = interp(v.env,m1[.!mask])
    else
        r[:] = interp(v.env,m1)
    end

    return r
end



"""
    vfun(id::Int,iy::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Array{Envelope},p::Param)

Calculate the period `it`, `iy`-state, discrete choice `id`-specific value function. Avoids interpolation in credit constrained region by using the analytic form of the value function (no need to interpolate expected value function when on the lower bound of assets.)
"""
function vfun(id::Int,iy::Int,it::Int,c1::Vector{Float64},m1::Vector{Float64},en::Array{Envelope,N} where N,p::Param)

    # L = en.L[id]
    v = en[id,iy,it]

    # computes v_{it}(m) = u(c) + beta v_{it+1}(m1)

    if length(getx(v)) < 2
        error("need more than 2 points in envelope object")
    end

    r = fill(NaN,size(m1))
    mask = m1.<getx(v)[2]
    mask = it==p.nT ? trues(size(mask)) : mask

    if all(mask)
        # in the credit constrained region:
        r[:] = u(c1,id==2,p) .+ p.beta * bound(v)
    elseif any(mask)
        r[mask] = u(c1[mask],id==2,p) .+ p.beta * bound(v)
        # elsewhere
        r[.!mask] = interp(v.env,m1[.!mask])
    else
        r[:] = interp(v.env,m1)
    end

    return r
end

"""
    ccp(m::Model,p::Param)

Conditional Choice probability of working
"""
function ccp(x::Matrix,p::Param) 
    #choice probability of the first row in 2-row matrix
    mx = maximum(x,dims = 1)
    mxx = x.-repeat(mx,size(x)[1],1)   # center values at max for numerical stability
    exp.(mxx[2,:]./p.lambda)./sum(exp.(mxx./p.lambda),dims = 1)[:]
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

Main body of the DC-EGM algorithm
"""
function dc_EGM!(m::Model,p::Param)
    for it in p.nT:-1:1
        println()
        @info("period = $it")

		if it==p.nT
            for id in 1:p.nD
    			# final period: consume everyting.
                # set the consumption function
                # remember this is a `MLine`, i.e. it has an x and a y Vector
                # x: endogenous grid m
                # y: optimal consumption at that grid x
                m.c[id,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,p.a_high)) )

                # initialize value function with vf(1) = 0
                m.v[id,it] = Envelope(MLine(vcat(p.a_lowT,p.a_high),vcat(0.0,NaN)) )
                # note that 0.0 as first value of the vfun is not innocuous here!
            end

		else
            for id in 1:p.nD
                working = id==2
                @info("id: $id")
                println()
    			# previous periods

    			# precomputed next period's cash on hand on all income state
    			# what's next period's cash on hand given you work/not today?
    			mm1 = m.m1[it+1][id]
                # println("mm1 = ")
                # display(mm1)

                # interpolate all d-choice next period's consumtion function on next cash on hand
                # println("next periods cons functions")
                # println([m.c[i,it+1].env for i in 1:2])
                c1 = interp([m.c[i,it+1].env for i in 1:2],mm1[:])
                c1[c1.<p.cfloor] .= p.cfloor   # no negative consumption
                # println("interpolated on next periods cash on hand. c1 = ")
                # display(c1)
                       
                # get next period's conditional value functions
                # as a matrix where each row is another discrete choice
                v1 = Matrix(hcat([vfun(jd,it+1,c1[jd,:],mm1[:],m.v,p) for jd in 1:p.nD]...)')

                # println("v1 = ")
                # display(v1[:,1:10])

                # get ccp to be a worker
                pwork = ccp(v1,p) * working

                # get marginal utility of that consumption
                mu1 = reshape(pwork .* up(c1[2,:],p) .+ (1 .- pwork) .* up(c1[1,:],p),p.na,p.ny)
                # println("mu1 = ")
                # display(mu1[1:10,:])

                # get expected marginal value of saving: RHS of euler equation
                # beta * R * E[ u'(c_{t+1}) ] 
                RHS = p.beta * p.R * mu1 * m.ywgt
                # println("RHS = $(RHS[1:10])")

                # optimal consumption today: invert the RHS of euler equation
                c0 = iup(RHS,p)

                # set optimal consumption function today. endo grid m and cons c0
                cline = MLine(m.avec .+ c0, c0)
                # store
                m.c[id,it] = Envelope(cline)
                # consumption function done.



    			# compute value function
    			# ----------------------

                if working
                    ev = reshape(logsum(v1,p),size(mm1)) * m.ywgt
                else
                    ev = reshape(vfun(1,it+1,c1[1,:],mm1[:],m.v,p),size(mm1)) * m.ywgt
                end
                vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev)

                # println(vline)


                # vline and cline may have backward-bending regions: let's prune those
                # SECONDARY ENVELOPE COMPUTATION

                if id==2   # only for workers
                    minx = minimum(vline.x)
                    if minx < vline.x[1]
                        # non-convex region lies inside credit constraint.
                        # endogenous x grid bends back before the first x grid point.
                        x0 = range(minx,stop = vline.x[1],length = floor(p.na/10)) # some points to the left of first x point
                        x0 = x0[1:end-1]
                        y0 = u(x0,working,p) + p.beta * ev[1]
                        prepend!(vline,x0,y0)
                        prepend!(cline,x0,x0)  # cons policy in credit constrained is 45 degree line
                    end

                    # split the vline at potential backward-bending points
                    # and save as Envelope object
                    m.v[id,it] = splitMLine(vline)  # splits line at backward bends

                    # if there is just one line (i.e. nothing was split in preceding step)
                    # then this IS a valid envelope
                    # else, need to compute the upper envelope.
                    if !m.v[id,it].env_set
                        # error()
                        upper_env!(m.v[id,it])   # compute upper envelope of this
                        # now need to clean the policy function as well
                        removed = getr(m.v[id,it])
                        for r in 1:length(removed)
                            if length(removed[r]) > 0
                                for ir in removed[r]
                                    delete!(m.c[id,it].L[r],ir.i)   # delete index i.r
                                end
                            end
                        end
                        # insert new intersections into consumption function
                        isecs = gets(m.v[id,it])
                        if length(isecs) > 0
                            for isec in 1:length(isecs)
                                I = isecs[isec]

                                # if that intersection is a new point
                                # i.e. intersection was not a member of any `MLine`
                                if I.new_point
                                    # insert intersection into env over cons function
                                    insert!(m.c[id,it].env,I.x,interp(m.c[id,it].env,[I.x]),I.i)

                                    # add to both adjacent `MLine` segments:
                                    # 1) append to end of segment preceding intersection:
                                    newy = interp(m.c[id,it].L[I.i],[I.x])
                                    append!(m.c[id,it].L[I.i],I.x,newy)
                                    # 1) prepend to beginning of segment following intersection:
                                    prepend!(m.c[id,it].L[I.i+1],I.x,newy)
                                end
                            end
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
                prepend!(m.c[id,it].env,p.a_low,0.0)
                # prepend!(m.v[id,it].env,p.a_low,ev[1])
                # do NOT prepend the value function with the special value from above.
    		end   # loop over discrete choice
    	end    # if final perio
    end     # loop over time
end


"""
    dcegm!(m::Model2,p::Param)

Main body of the DC-EGM algorithm version 2
"""
function dc_EGM!(m::Model2,p::Param)
    for it in p.nT:-1:1
        println()
        @info("period = $it")

        if it==p.nT
            for iy in 1:p.ny
                for id in 1:p.nD
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
                @info("iy: $iy")
                for id in 1:p.nD
                    working = id==2
                    @info("id: $id")
                    # previous periods

                    # precomputed next period's cash on hand on all income states
                    # what's next period's cash on hand given you work/not today?
                    mm1 = m.m1[it+1][id]
                    # println("mm1 = ")
                    # display(mm1)

                    # interpolate all d-choice next period's consumtion function on next cash on hand
                    # println("next periods cons functions")
                    # println([m.c[i,it+1].env for i in 1:2])
                    c1 = interp([m.c[i,iy,it+1].env for i in 1:2],mm1[:])
                    c1[c1.<p.cfloor] = p.cfloor   # no negative consumption
                    # println("interpolated on next periods cash on hand. c1 = ")
                    # display(c1)
                           
                    # get next period's conditional value functions
                    # as a matrix where each row is another discrete choice
                    v1 = Matrix(hcat([vfun(jd,iy,it+1,c1[jd,:],mm1[:],m.v,p) for jd in 1:p.nD]...)')

                    # println("v1 = ")
                    # display(v1[:,1:10])

                    # get ccp to be a worker
                    pwork = ccp(v1,p) * working


                    # get marginal utility of that consumption
                    mu1 = reshape(pwork .* up(c1[2,:],p) .+ (1 .- pwork) .* up(c1[1,:],p),p.na,p.ny)
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
                    cline = MLine(m.avec .+ c0, c0)
                    # store
                    m.c[id,iy,it] = Envelope(cline)
                    # consumption function done.



                    # compute value function
                    # ----------------------
                    if any(isnan.(logsum(v1,p)))
                        # @enter logsum(v1,p)
                    end

                    if working
                        ev = reshape(logsum(v1,p),size(mm1)) * m.ywgt[:,iy]
                    else
                        ev = reshape(vfun(1,iy,it+1,c1[1,:],mm1[:],m.v,p),size(mm1)) * m.ywgt[:,iy]
                    end
                    vline = MLine(m.avec .+ c0, u(c0,id==2,p) .+ p.beta * ev)

                    # println(vline)

                    if any(isnan.(ev)) 
                        println("ev = ")
                        display(ev)
                    end



                    # vline and cline may have backward-bending regions: let's prune those
                    # SECONDARY ENVELOPE COMPUTATION

                    if id==2   # only for workers
                        minx = minimum(vline.x)
                        if minx < vline.x[1]
                            # non-convex region lies inside credit constraint.
                            # endogenous x grid bends back before the first x grid point.
                            x0 = range(minx,stop = vline.x[1],length = floor(p.na/10)) # some points to the left of first x point
                            x0 = x0[1:end-1]
                            y0 = u(x0,working,p) + p.beta * ev[1]
                            prepend!(vline,x0,y0)
                            prepend!(cline,x0,x0)  # cons policy in credit constrained is 45 degree line
                        end

                        # split the vline at potential backward-bending points
                        # and save as Envelope object
                        m.v[id,iy,it] = splitMLine(vline)  # splits line at backward bends

                        # if there is just one line (i.e. nothing was split in preceding step)
                        # then this IS a valid envelope
                        # else, need to compute the upper envelope.
                        if !m.v[id,iy,it].env_set
                            upper_env!(m.v[id,iy,it])   # compute upper envelope of this
                            # println(m.c[id,iy,it].env)

                            # TODO
                            # bug
                            # this is not working properly. 
                            # for some reason i cannot get the set of removed points right
                            # here.
                            removed!(m.v[id,iy,it])
                            remove_c!(m.v[id,iy,it],m.c[id,iy,it])
                            # plot(m.c[id,iy,it].env)
                            # gui()

                            @assert(issorted(getx(m.v[id,iy,it])))
                            # display(hcat(getx(m.v[id,iy,it]),getx(m.c[id,iy,it])))
                            @assert(issorted(getx(m.c[id,iy,it])))
                            # insert new intersections into consumption function
                            isecs = gets(m.v[id,iy,it])
                            if length(isecs) > 0
                                for isec in 1:length(isecs)
                                    I = isecs[isec]

                                    # if that intersection is a new point
                                    # i.e. intersection was not a member of any `MLine`
                                    if I.new_point
                                        # insert intersection into env over cons function
                                        println("I.x = $(I.x)")
                                        println("I.i = $(I.i)")

                                        if !issorted(m.c[id,iy,it].env.x)
                                            println("m.c[id,iy,it].env.x = $(m.c[id,iy,it].env.x)")
                                            println("m.c[id,iy,it].env.y = $(m.c[id,iy,it].env.y)")
                                        end

                                        insert!(m.c[id,iy,it].env,I.x,interp(m.c[id,iy,it].env,[I.x]),I.i)

                                        # add to both adjacent `MLine` segments:
                                        # 1) append to end of segment preceding intersection:
                                        newy = interp(m.c[id,iy,it].L[I.i],[I.x])
                                        append!(m.c[id,iy,it].L[I.i],I.x,newy)
                                        # 1) prepend to beginning of segment following intersection:
                                        prepend!(m.c[id,iy,it].L[I.i+1],I.x,newy)
                                    end
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
                    prepend!(m.c[id,iy,it].env,p.a_low,0.0)
                    # prepend!(m.v[id,it].env,p.a_low,ev[1])
                    # do NOT prepend the value function with the special value from above.
                end # discrete choice
            end   # iy
        end    # if final perio
    end     # loop over time
end


function run2()
    p = Param()
    m = Model2(p)
    dc_EGM!(m,p)
    return (m,p)
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

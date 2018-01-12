

"""
    vfun(id::Int,it::Int,x::Vector{Float64},L::Line,p::Param)

Calculate the period `it`, discrete choice `id`-specific value function. Avoids interpolation in credit constrained region by using the analytic form of the value function (no need to interpolate expected value function when on the lower bound of assets.)
"""
function vfun(id::Int,x::Vector{Float64},en::Envelope,p::Param)

    # L = en.L[id]

    if length(en.L) < 2
        error("need more than 2 points in envelope object")
    end

    r = fill(NaN,size(x))
    mask = x.<en.L[id].x[2]
    mask = it==p.nT ? trues(mask) : mask

    if any(mask)
        # in the credit constrained region:
        r[mask] = u(x[mask],p,id==1) + p.beta * en.L[id].y[1]

        # elsewhere
        r[.!mask] = interp(en.L[id],x[.!mask])
    else
        r[:] = interp(en.L[id],x)
    end
end

"""
    ccp(m::Model,p::Param)

Conditional Choice probability of working
"""
function ccp(x::Matrix,p::Param) 
    #choice probability of the first row in multirow matrix
    mx = maximum(x,1)
    mxx = x.-repmat(mx,size(x)[1],1)   # center values at max for numerical stability
    exp.(mxx[1,:]./p.lambda)./sum(exp.(mxx./p.lambda),1)[:]
end

"""
    logsum(m::Model,p::Param)

Logsum of conditional values used in Expected value function.
"""
function logsum(x::Matrix,p::Param) 
    mx = maximum(x,1)
    mxx = x.-repmat(mx,size(x)[1],1)
    mx .+ p.lambda * log.( sum(exp.(mxx./p.lambda),1) )
end

"""
    dcegm!(m::Model,p::Param)

Main body of the DC-EGM algorithm
"""

function dc_EGM!(m::Model,p::Param)
    for it in p.nT:-1:1
    	for id in 1:p.nD
    		@debug(logger,"id: $id")
    		working = id==1

    		if it==p.nT
    			# final period: consume everyting.
                # set the consumption function
                # remember this is a `Line`, i.e. it has an x and a y Vector
                # x: endogenous grid m
                # y: optimal consumption at that grid x
                m.c[it][id] = Line(vcat(p.a_lowT,p.a_high),vcat(0.0,p.a_high))

                # initialize value function with vf(1) = 0
                m.v[it][id] = Line(vcat(p.a_lowT,p.a_high),vcat(0.0,NaN))

    		else
    			# previous periods

    			# precomputed next period's cash on hand on all income states
    			# what's next period's cash on hand given you work/not today?
    			mm1 = m.m1[it+1][id]

                # get next period's conditional value functions
                # as a matrix where each row is another discrete choice
                v1 = hcat([vfun(jd,mm1,m.v[it+1],p) for jd in 1:p.nD]...)'

                # get ccp to be a worker
                pwork = ccp(v1) * working

                # interpolate all d-choice next period's consumtion function on next cash on hand
                c1 = interp(m.c[it+1],mm1)
                c1[c1.<p.cfloor] = p.cfloor   # no negative consumption

                # get marginal utility of that consumption
                mu1 = pwork .* up(c1[1,:],p) .+ (1-pwork) .* up(c1[2,:],p)

                # get expected marginal value of saving: RHS of euler equation
                # beta * R * E[ u'(c_{t+1}) ] 
                RHS = p.beta * p.R * m.ywgt' * mu1

                # optimal consumption today: invert the RHS of euler equation
                c0 = iup(RHS,p)

                # set optimal consumption function today. endo grid m and cons c0
                m.c[it][id] = Line(m.avec .+ c0, c0)

    			# compute value function
    			# ----------------------

                if working
                    ev = m.ywgt' * repmat(logsum(v1,p),p.ny,1)
                else
                    ev = m.ywgt' * vfun(2,mm1,m.v[it+1],p)
                end
                m.v[it][id] = Line(m.avec .+ c0, u(c0,id,p) + p.beta * ev)
                # this value function may have backward-bending regions: let's prune those

                if id==1
                    minx = minimum(m.v[it][id].x)
                    if minx < m.v[it][id].x[1]
                        # non-convex region lies inside credit constraint.
                        # endogenous x grid bends back before the first x grid point.
                        x0 = linspace(minx,m.v[it][id].x[1],floor(p.na/10)) # some points to the left of first x point
                        x0 = x0[1:end-1]
                        y0 = u(x0,working,p) + p.beta * ev[1]
                        prepend!(m.v[it][id],x0,y0)
                    end
                    env = create_envelope(m.v[it][id])
                    upper_env!(env)   # compute upper envelope of this
                end

                # now need to clean the policy function as well
                if length(m.v[it].removed) > 0
                    # 1. remove m.v[it].removed from m.c[it]
                    # 2. i.e. from m.c[it].env and all m.c[it].L objects
                    for idel in 
                    # for isec in m.v[it].isects
                        # for each intersection in cond value functions
                        # drop the indices of removed points

                    end

                end

    		end  # if final period
    	end  # loop over discrete choice
    end
end
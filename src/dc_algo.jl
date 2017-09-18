

"""
    dcegm!(m::Model,p::Param)

Main body of the DC-EGM algorithm
"""

function d_EGM!(m::iidDModel,p::Param)
    for it in p.nT:-1:1
    	for id in 1:p.nD
    		@debug(logger,"id: $id")
    		working = convert(Bool,id - 1)

    		if it==p.nT
    			# final period: consume everyting.
                # set the consumption function
                set!(m.c[it],Line(vcat(0.0,p.a_high),vcat(0.0,p.a_high)))

                # set the value function
    			# this line assumes that you cannot die in debt
    			set_vbound!(m.v,it,id,0.0)
    			set_vbound!(m.v,it,0.0)
    			# v0 = u(c0,working,p) + p.beta * 0.0
    		else
    			# previous periods
    			# set values on lower bound
    			set_vbound!(m.c,it,id,0.0)
    			set_vbound!(m.m,it,id,0.0)
    			set_vbound!(m.v,it,id,NaN)

    			# precomputed next period's cash on hand on all income states
    			# what's next period's cash on hand given you work/not today?
    			mm1 = m.m1[it+1][id]

    			# next period's endogenous grid and cons function
    			# on the envelope!
    			# this must include the lower bound on both m and c
    			m1 = envvbound(m.m,it+1)  # p.na + 1
    			c1 = envvbound(m.c,it+1) 

    			# prepend next period's optimal policies with a zero to capture credit constraint
    			for ia in 1:p.na
    				for iy in 1:p.ny
    					tmpc = linearapprox(m1,c1,mm1[ia+p.na*(iy-1)])
    					m.c1[ia+p.na*(iy-1)] = max(tmpc,p.cfloor)
    				end
    			end

    			# get expected marginal value of saving: RHS of euler equation
    			# beta * R * E[ u'(c_{t+1}) ] 
    			Eu = up(m.c1,p) * m.ywgt
    			rhs = p.R * p.beta .* Eu

    			# set optimal consumption today from euler equation: invert marginal utility
    			set!(m.c,it,id,iup(Eu,p))
    			# set endogenous grid today
    			set!(m.m,it,id,cond(m.c,it,id) .+ m.avec)


    			# compute value function
    			# ======================

    			vv   = env(m.v,it+1)
    			tmpx = env(m.m,it+1)

    			# expected value function (na,ny)
    			fill!(m.ev,NaN)
    			# mask off values where we don't have to interpolate
    			if it==(p.nT-1)
    				mask = trues(size(mm1))
    			else
    				mask = mm1 .< env(m.m,it+1)[1]	# wherever potential next period's cash on hand (m.m1) is less than the second lowest grid point of the endogenous grid next period (lowest is 0), the agent will be credit constrained and will be saving zero 
    				# in that region, can use analytic form of v(work)
    			end

    			# again, compute EV directly for the credit constrained cases.
    			# remember that EV is the envelope over 2 conditional vfuns

    			for ia in 1:p.na
    				for iy in 1:p.ny
    					idx = ia+p.na*(iy-1)
    					# if credit constrained and next period is last: will choose retire and zero savings!
    					if it==(p.nT-1)  # retired next period
    						m.ev[idx] = u(max(mm1[idx],p.cfloor),false,p) + p.beta * get_vbound(m.v,it+1)
    					else
    						# if credit constrained and next period is not last: will choose work and zero savings!
    						if mask[idx]
    							m.ev[idx] = u(max(mm1[idx],p.cfloor),true,p) + p.beta * get_vbound(m.v,it+1)
    						else
    						# if not credit constrained: will save!
    							m.ev[idx] = linearapprox(tmpx,vv,mm1[idx])
    						end
    					end
    				end
    			end
    			# remember that mm1 is next period's cash on hand if TODAY discrete choice is id (e.g. work)
    			ev = m.ev * m.ywgt
    			# set value on bound
    			set_vbound!(m.v,it,id,ev[1])
    			# set remaining grid points of value function
    			set!(m.v,it,id, u(cond(m.c,it,id),working,p) + p.beta * ev)

    			# however, at this point there may be a problem in cond(m.v,it,id) because of secondary kinks coming in from env(m.v,it+1). Will now find wiggles where grid folds back onto itself.

    			if p.dorefinements

    				refine_grids!(m,it,id,p)
    				
    			end # if refinements
    		end  # if final period
    	end  # loop over discrete choice
    end
end
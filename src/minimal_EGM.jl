

function minimal_EGM(p::Param)
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    avec          = collect(range(0.0,p.a_high,length = p.na))
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    
    c[p.nT]       = [0.0,p.a_high]

    cg = cgrad(:viridis)
    cols = cg[range(0.0,stop=1.0,length=p.nT)]

    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=:topright,title="Consumption Function",
              xlims = (0,p.a_high),ylims = (0,p.a_high), color = cols[p.nT],
              xlab = "Cash on Hand", ylab = "Consumption")
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        Emu   = ywgt' * (1 ./ c1)   # Expected marginal utility (with log utility!). (p.na,1)
        rhs   = p.beta * p.R * Emu[:]   # RHS of euler equation
        c[it] = 1.0 ./ rhs   # inverse marginal utility function gives current consumption
        m[it] = avec .+ c[it]   # current period endogenous cash on hand grid. (p.na+1,1)

        # add credit constraint region
        # the lowest asset grid point on tomorrow's cash on hand (agrid[1]) corresponds to
        # saving exactly zero (or location on the lower bound of the asset grid)
        # well, we allow that people save zero (and consume all they have currently on hand).
        c[it] = vcat(0.0, c[it])   # prepend with 0
        m[it] = vcat(0.0, m[it])   #

        plot!(pl,m[it],c[it],label= it == 1 ? "$it" : "", color = cols[it])
    end
    pl = lens!(pl, [0, 2], [0, 2], inset = (1, bbox(0.2, 0.1, 0.25, 0.25)))
    pl
end

# artifial example: lowest m state corresponds to zero consumption, 
# even if that means negative cash on hand
# simple EGM that allows borrowing also in final period of life
# just shifts everything left
function minimal_EGM_borrow(p::Param)
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    avec          = collect(range(p.a_low,p.a_high,length = p.na))
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [p.a_low,p.a_high]    
    c[p.nT]       = [0.0,p.a_high]

    cg = cgrad(:viridis)
    cols = cg[range(0.0,stop=1.0,length=p.nT)]

    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=:topright,title="Consumption Function",
              xlims = (p.a_low,p.a_high),ylims = (0.0,p.a_high), color = cols[p.nT],
              xlab = "Cash on Hand", ylab = "Consumption")
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        Emu   = ywgt' * (1 ./ c1)   # Expected marginal utility (with log utility!). (p.na,1)
        rhs   = p.beta * p.R * Emu[:]   # RHS of euler equation
        c[it] = 1.0 ./ rhs   # inverse marginal utility function gives current consumption
        m[it] = avec .+ c[it]   # current period endogenous cash on hand grid. (p.na+1,1)

        # add credit constraint region
        # the lowest asset grid point on tomorrow's cash on hand (agrid[1]) corresponds to
        # saving exactly zero (or location on the lower bound of the asset grid)
        # well, we allow that people save zero (and consume all they have currently on hand).
        c[it] = vcat(0.0, c[it])   # prepend with 0
        m[it] = vcat(p.a_low, m[it])   #

        plot!(pl,m[it],c[it],label= it == 1 ? "$it" : "", color = cols[it])
    end
    pl = lens!(pl, [p.a_low, 2], [0.0, 2], inset = (1, bbox(0.2, 0.1, 0.25, 0.25)))
    m,c,pl
end

# artifial example: lowest m state corresponds to zero consumption, 
# even if that means negative cash on hand
# simple EGM that allows borrowing also in final period of life
# just shifts everything left
# function minimal_EGM_borrow(p::Param)
#     nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
#     yvec          = sqrt(2.0) * p.sigma .* nodes
#     ywgt          = weights .* pi^(-0.5)
#     avec          = collect(range(p.a_low,p.a_high,length = p.na))
#     m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
#     c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
#     m[p.nT]       = [0.0,p.a_high]    
#     c[p.nT]       = [0.0,p.a_high]

#     cg = cgrad(:viridis)
#     cols = cg[range(0.0,stop=1.0,length=p.nT)]

#     pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=:topright,title="Consumption Function",
#               xlims = (p.a_low,p.a_high),ylims = (0.0,p.a_high), color = cols[p.nT],
#               xlab = "Cash on Hand", ylab = "Consumption")
#     # cycle back in time
#     for it in p.nT-1:-1:1
#         w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
#         # get next period consumption on that wealth w1
#         # interpolate on next period's endogenous grid m[it+1].
#         # notice that the `interpolate` object needs to be able to extrapolate
#         c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
#         c1[c1.<0] .= p.cfloor     # don't allow negative consumption
#         Emu   = ywgt' * (1 ./ c1)   # Expected marginal utility (with log utility!). (p.na,1)
#         rhs   = p.beta * p.R * Emu[:]   # RHS of euler equation
#         c[it] = 1.0 ./ rhs   # inverse marginal utility function gives current consumption
#         m[it] = avec .+ c[it]   # current period endogenous cash on hand grid. (p.na+1,1)

#         # add credit constraint region
#         # the lowest asset grid point on tomorrow's cash on hand (agrid[1]) corresponds to
#         # saving exactly zero (or location on the lower bound of the asset grid)
#         # well, we allow that people save zero (and consume all they have currently on hand).
#         c[it] = vcat(0.0, c[it])   # prepend with 0
#         m[it] = vcat(p.a_low, m[it])   #

#         plot!(pl,m[it],c[it],label= it == 1 ? "$it" : "", color = cols[it])
#     end
#     pl = lens!(pl, [p.a_low, 2], [0.0, 2], inset = (1, bbox(0.2, 0.1, 0.25, 0.25)))
#     m,c,pl
# end


function minimal_EGM_bequest(p::Param)
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))
    avec          = scaleGrid(0.0,p.a_high,p.na,logorder = 1)
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    v             = Vector{Float64}[Float64[] for i in 1:p.nT]   # value function on m

    m[p.nT]       = avec    # no debt in last period possible
    c[p.nT]       = avec    # consume all you have left
    v[p.nT]       = p.ν > 0 ? bequest(avec,p) : zeros(p.na)

    # plot setup
    cg = cgrad(:plasma)
    cols = cg[range(0.0,stop=1.0,length=p.nT)]
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,title="Consumption",color = cols[p.nT], xlims = (0,p.a_high),ylims = (0,p.a_high))
    pv = plot(m[p.nT],v[p.nT],label="$(p.nT)",leg=:bottomright,title="Values",color = cols[p.nT], xlims = (0,p.a_high), ylims = (-20,2))

    for it in (p.nT-1):-1:1
        w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption

        # take care of bequest: if this is the penultimate period,
        # need to decide how much to leave for bequest
        if it == p.nT-1
            if p.ν > 0
                # if there is a bequest motive at all
                # this is like any other savings decision, except for that
                # next (i.e. final) period utility is the bequest function.
                # this looks like u(b + bbar) where bbar > 0 governs says that people can leave
                # zero bequests and survive (i.e. do not get u(0) = -Inf).
                mu1 = ywgt' * (1 ./ (c1 .+ p.bbar))  # expected marginal utility of consumption next period
                rhs = p.beta * p.R * p.ν * mu1[:]   # RHS of Euler equation
                c[it] = vcat(0.0, (1.0 ./ rhs[:])...)   # current period consumption vector. (p.na+1,1)
            else
                # if there is no bequest motive, this is like a standard period
                mu1 = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
                rhs = p.beta * p.R * mu1[:]
                c[it] = vcat(0.0, (1.0 ./ rhs[:])...)   # current period consumption vector. (p.na+1,1)
            end
        else
            mu1 = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
            rhs = p.beta * p.R * mu1[:]
            c[it] = vcat(0.0, (1.0 ./ rhs[:])...)   # current period consumption vector. (p.na+1,1)
        end

        # we prepend c and m with (0) here to create the borrowing constrained region
        m[it] = vcat(0.0, avec .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        plot!(pl,m[it],c[it],label = it == 1 ? "$it" : "",color = cols[it])

        # compute value function
        # next period values at all states w1
        v1 = reshape(extrapolate(interpolate((m[it+1],),v[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        ev = ywgt' * v1  # integrate out wage shocks
        v[it] = vcat(ev[1], u(c[it][2:end],p) .+ p.beta * ev[:])  # ev[1] is expected value of saving exactly zero.
        plot!(pv, m[it][2:end], v[it][2:end],label = it == 1 ? "$it" : "",color = cols[it])

    end

    pl = lens!(pl, [0, 2], [0, 2], inset = (1, bbox(0.2, 0.05, 0.3, 0.25)))
    return plot(pl,pv,layout = (1,2))
end




function minimal_EGM()
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))
    avec          = scaleGrid(0.0,p.a_high,p.na,logorder = 1)
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,title="0 borrowing")
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
        c[it] = vcat(0.0, 1.0 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector. (p.na+1,1)
        m[it] = vcat(0.0, avec .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        plot!(pl,m[it],c[it],label="$it")
    end
    return (m,c,pl)
end

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
    v[p.nT]       = bequest(avec,p)

    # plot setup
    cg = cgrad(:plasma)
    cols = cg[range(0.0,stop=1.0,length=p.nT)]
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,title="Consumption",color = cols[p.nT], xlims = (0,p.a_high),ylims = (0,p.a_high))
    pv = plot(m[p.nT],v[p.nT],label="$(p.nT)",leg=:bottomright,title="Values",color = cols[p.nT], xlims = (0,p.a_high), ylims = (-20,-1))
    # cycle back in time
    for it in p.nT-1:-1:1
        w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec'   # w1 = y + yshock*R*savings:  next period wealth at all states. (p.ny,p.na)
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)

        # we prepend c and m with (0) here to create the borrowing constrained region
        c[it] = vcat(0.0, 1.0 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector. (p.na+1,1)
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

function minimal_EGM_neg(;alow = -1.0,brate = false)
    p             = Param()
    nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
    yvec          = sqrt(2.0) * p.sigma .* nodes
    ywgt          = weights .* pi^(-0.5)
    # avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))

    # η = NBL(yvec[1],p)  # compute natural borrowing limit: maximal borrowing if c=0 in all future periods and worst income
    # if blim
    #     avec          = [scaleGrid(η[it],p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
    #     push!(avec, scaleGrid(0.0,p.a_high,p.na,logorder = 1))  # last period
    # else
        avec          = [scaleGrid(alow,p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
        push!(avec, scaleGrid(0.0,p.a_high,p.na,logorder = 1))  # last period
    # end
    m             = Vector{Float64}[Float64[] for i in 1:p.nT]   # endogenous grid
    c             = Vector{Float64}[Float64[] for i in 1:p.nT]   # consumption function on m
    m[p.nT]       = [0.0,p.a_high]    # no debt in last period possible
    c[p.nT]       = [0.0,p.a_high]
    # m[p.nT]       = [p.a_low,0.0,p.a_high]    # no debt in last period possible
    # c[p.nT]       = [p.cfloor,p.cfloor,p.a_high]
    # ti = blim ? "Ntl. br:$brate" : "Fixed. br:$brate"
    pl = plot(m[p.nT],c[p.nT],label="$(p.nT)",leg=false,title = ti)
    # cycle back in time
    for it in p.nT-1:-1:1
        if !brate
            w1 = 0.0 .+ exp.(yvec) .+ p.R.*avec[it+1]'   # w1 = y + yshock + R*savings:  next period wealth at all states. (p.ny,p.na)
        else

            w1 = 0.0 .+ exp.(yvec) .+ vcat(avec[it+1][avec[it+1] .< 0 ].*(p.R*1.3),avec[it+1][avec[it+1] .>= 0 ].*(p.R))'  # w1 = y + yshock + R*savings:  next period wealth at all states. (p.ny,p.na)
        end
        # get next period consumption on that wealth w1
        # interpolate on next period's endogenous grid m[it+1].
        # notice that the `interpolate` object needs to be able to extrapolate
        c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Line())(w1[:]) ,p.ny,p.na)
        # if it == p.nT-3
            # println("mt+1 = $(m[it+1])")
            # println("wt+1 = $w1")
            # println("ct+1 = $c1")
        # end
        c1[c1.<0] .= p.cfloor     # don't allow negative consumption
        rhs = ywgt' * (1 ./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
        c[it] = vcat(0.0, 1 ./ (p.beta * p.R * rhs[:])...)   # current period consumption vector ∈ [0,∞]
        m[it] = vcat(avec[it][1], avec[it] .+ c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
        # plot!(pl,m[it],c[it],label="$it")
        plot!(pl,m[it],c[it])
    end
    return (m,c,pl)
end

function plot_minimal()
    p1 = minimal_EGM()
    p2 = minimal_EGM_neg()
    p3 = minimal_EGM_neg(brate = true)
    p4 = minimal_EGM_neg(blim = false)
    p5 = minimal_EGM_neg(blim = false,brate = true)
    plot(p1[3],p2[3],p3[3],p4[3],p5[3],layout = (1,5),xlims = (-10,10),size = (800,400))
end

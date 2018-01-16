
"""
Holds the user-set parameter values. 

**values set in param.json file:**

* `gamma`: CRRA coefficient
* `beta`: discount factor
* `Tbar`: maximal age
* `R`: gross interest rate (i.e. 1+r )
* `na`: number of grid points for assets
* `ny`: number of grid points for income
* `a_lowT`: lower bound on assets in final period
* `a_low`: lower bound on assets
* `a_high`: upper bound on assets
* `mu`: unconditional mean of income (iid case)
* `sigma`: unconditional variance of income (iid case)
* `rho_z`: AR1 coefficient of income (AR1 case)
* `eps_z`: standard deviation of AR1 innovation (AR1 case)
* `dorefinements`: boolean switch of whether to filter out whiggles
* `alpha`: disutility of work
* `inc0`: income equation
* `inc1`: income equation
* `inc2`: income equation
* `cfloor`: consumption floor


"""
mutable struct Param

	gamma                 :: Float64  # CRRA
	neg_gamma             :: Float64
	oneminusgamma         :: Float64
	oneover_oneminusgamma :: Float64
	neg_oneover_gamma     :: Float64
	beta                  :: Float64
	sigma                 :: Float64
	lambda                :: Float64
	R                     :: Float64
	na                    :: Int 	# asset grid
	ny                    :: Int   # income grid (support points)
	nT                    :: Int   # maximal age
	a_high                :: Float64
	a_low                 :: Float64
	a_lowT                :: Float64
	nD                    :: Int # number of discrete choices
	cfloor                :: Float64
	alpha                 :: Float64
	inc0                 :: Float64
	inc1                 :: Float64
	inc2                 :: Float64

	# constructor 
    function Param(;par::Dict=Dict())
		f=open(joinpath(dirname(@__FILE__),"param.json")) 
		j = JSON.parse(f)
		close(f)
    	this = new()
    	for (k,v) in j
    		setfield!(this,Symbol(k),v["value"])
    	end

        if length(par) > 0
            # override parameters from dict par
            for (k,v) in par
                setfield!(this,k,v)
            end
        end

    	# derived
		this.neg_gamma             = (-1.0) * this.gamma
		this.oneminusgamma         = 1.0 - this.gamma
		this.oneover_oneminusgamma = 1.0 / this.oneminusgamma
		this.neg_oneover_gamma     = (-1.0) / this.gamma

		return this
	end
end


"""
Binary Choice Model with iid income uncertainty

uses cash-on-hand m=y+a as state variable
"""
mutable struct Model

	# nD is number of discrete choices: nD = 2

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income weights

	# intermediate objects (na,ny,nD)
	m1::Dict{Int,Dict}	# a dict[it] for each period
	c1::Matrix{Float64}
	ev::Matrix{Float64}

	# result objects
	v :: Matrix{Envelope}  # a vector of Envelope objects
	c :: Matrix{Envelope}

	"""
	Constructor for discrete choice Model
	"""
	function Model(p::Param)

		this = new()
		# avec          = scaleGrid(0.0,p.a_high,p.na,2)
		this.avec          = collect(linspace(p.a_low,p.a_high,p.na))

		# fedors vesion:
		# nodes,weights = quadpoints(p.ny,0,1) 
		# N = Normal(0,1)
		# nodes = quantile.(N,nodes)
		# this.yvec = nodes * p.sigma
		# this.ywgt = weights

		# my version:
		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
		nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
		this.yvec = sqrt(2.0) * p.sigma .* nodes
		this.ywgt = weights .* pi^(-0.5)

		# precompute next period's cash on hand.
		# (na,ny,nD)
		# iD = 1: no work
		# iD = 2: work
		this.m1 = Dict(it => Dict(id => Float64[this.avec[ia]*p.R + income(it,p,this.yvec[iy]) * (id-1) for ia in 1:p.na, iy in 1:p.ny  ] for id=1:p.nD) for it=1:p.nT)
		this.c1 = zeros(p.na,p.ny)
		this.ev = zeros(p.na,p.ny)

		this.v = [Envelope(Line(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, it in 1:p.nT]
		this.c = [Envelope(Line(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, it in 1:p.nT]
		# dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

		return this
	end
end


function minimal_EGM(;dplot=false)
	p             = Param()
	nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
	yvec          = sqrt(2.0) * p.sigma .* nodes
	ywgt          = weights .* pi^(-0.5)
	avec          = collect(linspace(p.a_low,p.a_high,p.na))
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
		c1 = reshape(extrapolate(interpolate((m[it+1],),c[it+1],Gridded(Linear())),Linear())[w1[:]] ,p.ny,p.na)  
		c1[c1.<0] = p.cfloor
		rhs = ywgt' * (1./ c1)   # rhs of euler equation (with log utility!). (p.na,1)
		c[it] = vcat(0, 1./(p.beta * p.R * rhs[:])...)   # current period consumption vector. (p.na+1,1)
		m[it] = vcat(p.a_low, avec.+c[it][2:end]...)   # current period endogenous cash on hand grid. (p.na+1,1)
		if dplot
			plot!(m[it],c[it],label="$it")
		end
	end
	if dplot
		gui()
	end

	return (m,c)
end
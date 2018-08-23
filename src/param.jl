
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

	



"""
Model 2
"""
mutable struct Model2

	# nD is number of discrete choices: nD = 2

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # income weights

	# intermediate objects (na,ny,nD)
	m1::Dict{Int,Dict}	# a dict[it] for each period
	c1::Matrix{Float64}
	ev::Matrix{Float64}

	# result objects
	v :: Array{Envelope}  # arrays of Envelope objects
	c :: Array{Envelope}

	"""
	Constructor for discrete choice Model
	"""
	function Model2(p::Param)

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
		ywgt = weights .* pi^(-0.5)
		ywgt = ywgt .+ ywgt'
		this.ywgt = ywgt ./ sum(ywgt,2)


		# precompute next period's cash on hand.
		# (na,ny,nD)
		# iD = 1: no work
		# iD = 2: work
		this.m1 = Dict(it => Dict(id => Float64[this.avec[ia]*p.R + income(it,p,this.yvec[iy]) * (id-1) for ia in 1:p.na, iy in 1:p.ny  ] for id=1:p.nD) for it=1:p.nT)
		this.c1 = zeros(p.na,p.ny)
		this.ev = zeros(p.na,p.ny)

		this.v = [Envelope(Line(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, iy in 1:p.ny, it in 1:p.nT]
		this.c = [Envelope(Line(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, iy in 1:p.ny,it in 1:p.nT]
		# dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

		return this
	end
end

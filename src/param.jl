
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
struct Param

	gamma                 :: Float64  # CRRA
	neg_gamma             :: Float64
	oneminusgamma         :: Float64
	oneover_oneminusgamma :: Float64
	neg_oneover_gamma     :: Float64
	beta                  :: Float64
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
		this.oneover_oneminusgamma = 1.0 / oneminusgamma
		this.neg_oneover_gamma     = (-1.0) / this.gamma

		return this
	end
end


"""
Binary Choice Model with iid income uncertainty

uses cash-on-hand m=y+a as state variable
"""
mutable struct iidDModel

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
	m :: Dict{Int,Envelope}
	v :: Dict{Int,Envelope}
	c :: Dict{Int,Envelope}

	"""
	Constructor for iid Dchoice Model
	"""
	function iidDModel(p::Param)

		this = new()
		# avec          = scaleGrid(0.0,p.a_high,p.na,2)
		this.avec          = collect(linspace(p.a_low,p.a_high,p.na))
		# nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
		nodes,weights = quadpoints(p.ny,0,1)  # from FastGaussQuadrature
		N = Normal(0,1)
		nodes = quantile(N,nodes)

		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
		# yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		this.yvec = nodes * p.sigma
		# ywgt = weights .* pi^(-0.5)
		this.ywgt = weights

		# precompute next period's cash on hand.
		# (na,ny,nD)
		# iD = 1: no work
		# iD = 2: work
		this.m1 = [it => [id => Float64[avec[ia]*p.R + income(it,p,yvec[iy]) * (id-1) for ia in 1:p.na, iy in 1:p.ny  ] for id=1:p.nD] for it=1:p.nT]
		this.c1 = zeros(p.na,p.ny)
		this.ev = zeros(p.na,p.ny)

		# dicts
		# m = [it => Envelope([id => zeros(p.na) for id in 1:2],[id => 0.0 for id in 1:2], 0.0, zeros(p.na)) for it in 1:p.nT]
		this.m = [it => Envelope(0.0) for it in 1:p.nT]
		this.v = [it => Envelope(0.0) for it in 1:p.nT]
		this.c = [it => Envelope(0.0) for it in 1:p.nT]
		# dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

		return this
	end
end
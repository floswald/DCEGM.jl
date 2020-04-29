
abstract type Model end




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
	cfloor_plot           :: Float64
	alpha                 :: Float64
	alphaT                 :: Float64
	inc0                 :: Float64
	inc1                 :: Float64
	inc2                 :: Float64
	ρ                    :: Float64
	k                    :: Int

	# simulation
	nsims                 :: Int64
	initw0                :: Float64   # low/high bounds on initial wealth from this interval
	initw1                :: Float64   # low/high bounds on initial wealth from this interval
	rseed                 :: Int64     # random generator seed

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


		# this.beta = 1/this.R  # hard wire that beta = 1/R
		# this is not good.

		return this
	end
end


# """
# Binary Choice GModel with iid income uncertainty

# uses cash-on-hand m=y+a as state variable
# """
# mutable struct GModel

# 	# nD is number of discrete choices: nD = 2

# 	# computation grids
# 	avec::Vector{Float64}
# 	yvec::Vector{Float64}   # income support
# 	ywgt::Vector{Float64}   # income weights

# 	# intermediate objects (na,ny,nD)
# 	m1::Dict{Int,Dict}	# a dict[it] for each period
# 	c1::Matrix{Float64}
# 	ev::Matrix{Float64}

# 	# result objects
# 	v :: Matrix{Envelope}  # a vector of Envelope objects
# 	c :: Matrix{Envelope}

# 	"""
# 	Constructor for discrete choice GModel
# 	"""
# 	function GModel(p::Param)

# 		this = new()
# 		# avec          = scaleGrid(0.0,p.a_high,p.na,2)
# 		this.avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))

# 		# fedors vesion:
# 		# nodes,weights = quadpoints(p.ny,0,1)
# 		# N = Normal(0,1)
# 		# nodes = quantile.(N,nodes)
# 		# this.yvec = nodes * p.sigma
# 		# this.ywgt = weights

# 		# my version:
# 		# for y ~ N(mu,sigma), integrate y with
# 		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
# 		nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
# 		this.yvec = sqrt(2.0) * p.sigma .* nodes
# 		this.ywgt = weights .* pi^(-0.5)

# 		# precompute next period's cash on hand.
# 		# (na,ny,nD)
# 		# iD = 1: no work
# 		# iD = 2: work
# 		this.m1 = Dict(it => Dict(id => Float64[this.avec[ia]*p.R + income(it,p,this.yvec[iy]) * (id-1) for ia in 1:p.na, iy in 1:p.ny  ] for id=1:p.nD) for it=1:p.nT)
# 		this.c1 = zeros(p.na,p.ny)
# 		this.ev = zeros(p.na,p.ny)

# 		this.v = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, it in 1:p.nT]
# 		this.c = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, it in 1:p.nT]
# 		# dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

# 		return this
# 	end
# end



"""
Fedor's Model
"""
mutable struct FModel <: Model
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income support

	# intermediate objects (na,nD)
	m1::Dict{Int,Dict}	# a dict[it] for each period
	c1::Matrix{Float64}
	ev::Matrix{Float64}

	# result objects
	v :: Array{Envelope}  # arrays of Envelope objects
	c :: Array{Envelope}


	function FModel(p::Param)

		this = new()

		# fedors version:
		nodes,weights = quadpoints(p.ny,0,1)
		N = Normal(0,1)
		nodes = quantile.(N,nodes)
		this.yvec = nodes * p.sigma
		this.ywgt = weights

		# simulation shocks

		this.avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))

		# precompute next period's cash on hand.
		# (na,ny,nD)
		# iD = 1: tomorrow work
		# iD = 2: tomorrow no work - absorbing state and retire
		# notice: you decide to work today (t), but your shock realises tomorrow (t+1) in their formulation
		this.m1 = Dict(it => Dict(id => Float64[this.avec[ia]*p.R .+ income(it+1,p,this.yvec[iy]) * (id==1) for iy in 1:p.ny , ia in 1:p.na ] for id=1:p.nD) for it=1:(p.nT-1))

		# result arrays: matrices of type Envelope.
		this.v = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, it in 1:p.nT]
		this.c = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD ,it in 1:p.nT]

		return this
	end
end

"""
Fedor's Model with bankruptcy
"""
mutable struct BModel <: Model
	avec::Vector{Float64}
	aposvec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # income support

	# intermediate objects (na,nD)
	m1::Dict{Int,Dict}	# a dict[it] for each period
	# result objects
	v :: Array{Envelope}  # arrays of Envelope objects
	c :: Array{Envelope}
	vbk :: Array{Envelope}  # arrays of Envelope objects
	cbk :: Array{Envelope}
	iazero :: Int  # index of first non-negative asset state



	function BModel(p::Param)

		this = new()

		if p.ρ == 0
			nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
			this.yvec = sqrt(2.0) * p.sigma .* nodes
			this.ywgt = reshape(repeat(weights .* pi^(-0.5),inner = p.ny),p.ny,p.ny)  # make a matrix
		else
			# version with income persistence
			this.yvec, this.ywgt = rouwenhorst(p.ρ,0,p.sigma,p.ny)
		end

		if p.a_low >= 0
			error("need a_low < 0 for bankruptcy model")
		end

		# simulation shocks

		this.avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))
		this.aposvec       = collect(range(0.0,stop = p.a_high,length = p.na))
		this.iazero = findfirst(this.avec .>= 0)

		# precompute next period's cash on hand.
		# (na,ny,nState)
		# state = 1: tomorrow bk flag off
		# state = 2: tomorrow bk flag on
		this.m1 = Dict(it => Dict(id =>
		                          Float64[p.R* (this.aposvec[ia]*(id==2) + (1 - (id==2))*this.avec[ia]).+ income(it,p,this.yvec[iy]) for iy in 1:p.ny , ia in 1:p.na] for id=1:(p.nD)) for it=2:(p.nT))

		# result arrays: matrices of type Envelope.
		# this allocation is only to reserve about the right amount of memory. those will be overwritten in the algo.
		this.v = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, iy in 1:p.ny, it in 1:p.nT]
		this.c = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, iy in 1:p.ny,it in 1:p.nT]
		# for bankruptcy flag on there is no discrete choice
		this.vbk = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for iy in 1:p.ny, it in 1:p.nT]
		this.cbk = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for iy in 1:p.ny,it in 1:p.nT]


		return this
	end
end



"""
General GModel
"""
mutable struct GModel <: Model

	# nD is number of discrete choices: nD = 2

	# computation grids
	avec::Vector{Float64}  # each period has its own avec
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # income weights

	# intermediate objects (na,ny,nD)
	m1::Dict{Int,Dict}	# a dict[it] for each period
	c1::Matrix{Float64}
	ev::Matrix{Float64}

	# result objects
	v :: Array{Envelope}  # arrays of Envelope objects
	c :: Array{Envelope}

	# simulation settings
	wshocks :: Matrix{Float64}  # wage shocks
	dshocks :: Matrix{Float64}  # discrete choice shocks
	w0shocks :: Vector{Float64}  # initial wealth shock
	y0shocks :: Vector{Int64}  # initial income state

	"""
	Constructor for discrete choice GModel
	"""
	function GModel(p::Param)

		this = new()

		# fedors version:
		# nodes,weights = quadpoints(p.ny,0,1)
		# N = Normal(0,1)
		# nodes = quantile.(N,nodes)
		# this.yvec = nodes * p.sigma
		# # this.ywgt = weights'
		# this.ywgt = reshape(repeat(weights,inner = p.ny),p.ny,p.ny)  # make a matrix


		# my version:
		# for y ~ N(mu,sigma), integrate y with
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature

		if p.ρ == 0
			nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
			this.yvec = sqrt(2.0) * p.sigma .* nodes
			this.ywgt = reshape(repeat(weights .* pi^(-0.5),inner = p.ny),p.ny,p.ny)  # make a matrix
		else
			# version with income persistence
			this.yvec, this.ywgt = rouwenhorst(p.ρ,0,p.sigma,p.ny)
		end


		this.avec          = scaleGrid(p.a_low,p.a_high,p.na,logorder = 2)
		# this.avec          = [collect(range(p.a_low,stop = p.a_high,length = p.na))]
		# this.avec          = [scaleGrid(p.a_low,p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
		# this.avec          = [scaleGrid(η[it],p.a_high,p.na,logorder = 1) for it in 1:p.nT-1]
		# this.avec          = [scaleGrid(p.a_low,p.a_high,p.na,logorder = 0) for it in 1:p.nT-1]
		# push!(this.avec, scaleGrid(0.0,p.a_high,p.na,logorder = 0))  # last period
		# this.avec          = collect(range(p.a_low,stop = p.a_high,length = p.na))

		# precompute next period's cash on hand.
		# (na,ny,nD)
		# iD = 1: tomorrow work
		# iD = 2: tomorrow no work
		# this.m1 = Dict(it => Dict(id => Float64[this.avec[ia]*p.R .+ income(it,p,this.yvec[iy]) * (id==1) *(it < p.nT) for iy in 1:p.ny , ia in 1:p.na ] for id=1:p.nD) for it=2:p.nT)
		this.m1 = Dict(it => Dict(id => Float64[this.avec[ia]*((this.avec[ia] < 0)*0.5 + p.R)  .+ income(it,p,this.yvec[iy]) * (id==1) for iy in 1:p.ny , ia in 1:p.na ] for id=1:p.nD) for it=2:p.nT)
		this.c1 = zeros(p.na,p.ny)
		this.ev = zeros(p.na,p.ny)

		this.v = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, iy in 1:p.ny, it in 1:p.nT]
		this.c = [Envelope(MLine(fill(NaN,(p.na)),fill(NaN,(p.na)))) for id in 1:p.nD, iy in 1:p.ny,it in 1:p.nT]
		# dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

		return this
	end
end


function gmodel()
	p = Param()
	(GModel(p),p)
end

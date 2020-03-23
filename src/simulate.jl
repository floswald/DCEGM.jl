

mutable struct Simulation
	w0        :: Matrix{Float64}
	w1        :: Matrix{Float64}
	cons      :: Matrix{Float64}
	shocks    :: Matrix{Float64}
	inc       :: Matrix{Float64}
	prob_work :: Matrix{Float64}
	ret_age   :: Vector{Int64}
	worker    :: Matrix{Int64}
	function Simulation(p::Param)
		this = new()
		this.w0        = fill(NaN,p.nsims,p.nT)  # beginning of period wealth
		this.w1        = fill(NaN,p.nsims,p.nT)  # end of period wealth
		this.cons      = fill(NaN,p.nsims,p.nT)
		this.shocks    = fill(NaN,p.nsims,p.nT)
		this.inc       = fill(NaN,p.nsims,p.nT)
		this.worker    = zeros(Int,p.nsims,p.nT)
		this.prob_work = fill(NaN,p.nsims,p.nT)
		this.ret_age   = zeros(Int,p.nsims)
		this
	end
end

function sim(m::FModel,p::Param)

	s = Simulation(p)

	working   = trues(p.nsims)  # indicator of currentlcy working or not
	vmat      = zeros(2,p.nsims)

	N = Normal(0,1)  # normal dist N(0,σ)


	for it in 1:p.nT
		if it == 1
			# draw from initial distributions
			s.w0[:, it] = p.initw0 .+ rand(p.nsims)*(p.initw1 - p.initw0)
			s.worker[:, it] .= 1   # all work in frist period
		else

			# working state
			s.worker[ working, it] .= 1 # chose working, so stay worker
			s.worker[.!(working), it] .= 2 # chose retire
			s.ret_age[(s.worker[:,it-1] .== 1) .& (s.worker[:,it] .== 2)] .= it

			s.shocks[working,it] = quantile.(N,rand(sum(working))) .* p.sigma

			# income
			s.inc[.!(working), it] .= 0.0
			s.inc[ working, it] .= income(it,p,s.shocks[working,it])


			# end of period wealth in period it
			s.w0[ :     , it] .= s.w1[:      , it-1]*p.R
			s.w0[working, it] .= s.w1[working, it-1]*p.R .+ s.inc[ working, it]
		end
		# consumption
		s.cons[ working,it]    = gety(interp(m.c[1,it].env,s.w0[ working, it]))
		s.cons[.!(working),it] = gety(interp(m.c[2,it].env,s.w0[.!(working), it]))
		# end of period wealth
		s.w1[:,it] = s.w0[:,it] - s.cons[:,it]

		# CCP to remain worker
		vmat[1,:] = vfun(1,it, s.cons[:,it] ,s.w0[:,it], m.v[1,it],p )
		vmat[2,:] = vfun(2,it, s.cons[:,it] ,s.w0[:,it], m.v[2,it],p )
		s.prob_work[:,it] = (s.worker[:,it] .== 1) .* ccp(vmat,p)  # prob to remain worker

		# discrete choice
		working[:] = s.prob_work[:,it] .> rand(p.nsims)
	end
	s
end

function rs(;par = Dict())
	m,p = DCEGM.runf(par = par)
	sim(m,p)
end

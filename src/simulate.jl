

mutable struct Simulation
	w0        :: Matrix{Float64}
	w1        :: Matrix{Float64}
	cons      :: Matrix{Float64}
	ystate    :: Matrix{Int64}
	inc       :: Matrix{Float64}
	prob_work :: Matrix{Float64}
	ret_age   :: Vector{Int64}
	worker    :: Matrix{Int64}
	p :: Param
	function Simulation(p::Param)
		this = new()
		this.w0        = fill(NaN,p.nsims,p.nT)  # beginning of period wealth
		this.w1        = fill(NaN,p.nsims,p.nT)  # end of period wealth
		this.cons      = fill(NaN,p.nsims,p.nT)
		this.ystate    = zeros(Int,p.nsims,p.nT)
		this.inc       = fill(NaN,p.nsims,p.nT)
		this.worker    = zeros(Int,p.nsims,p.nT)
		this.prob_work = fill(NaN,p.nsims,p.nT)
		this.ret_age   = zeros(Int,p.nsims)
		this.p = deepcopy(p)
		this
	end
end

mutable struct BSimulation
	ainit        :: Matrix{Float64}
	w0        :: Matrix{Float64}
	w1        :: Matrix{Float64}
	cons      :: Matrix{Float64}
	ystate    :: Matrix{Int64}
	inc       :: Matrix{Float64}
	prob_nofile :: Matrix{Float64}
	file_age   :: Vector{Int64}
	bkflag    :: Matrix{Bool}  # on/off
	vfile :: Matrix{Float64}
	vnofile :: Matrix{Float64}
	p :: Param
	function BSimulation(p::Param)
		this = new()
		this.ainit        = fill(NaN,p.nsims,p.nT)  # beginning of period wealth
		this.w0        = fill(NaN,p.nsims,p.nT)  # beginning of period wealth
		this.w1        = fill(NaN,p.nsims,p.nT)  # end of period wealth
		this.cons      = fill(NaN,p.nsims,p.nT)
		this.vfile      = fill(NaN,p.nsims,p.nT)
		this.vnofile      = fill(NaN,p.nsims,p.nT)
		this.ystate    = zeros(Int,p.nsims,p.nT)
		this.inc       = fill(NaN,p.nsims,p.nT)
		this.bkflag    = falses(p.nsims,p.nT)
		this.prob_nofile = fill(NaN,p.nsims,p.nT)
		this.file_age   = zeros(Int,p.nsims)
		this.p = deepcopy(p)
		this
	end
end

function sim(m::FModel,p::Param)

	s = Simulation(p)

	working   = trues(p.nsims)  # indicator of currentlcy working or not
	vmat      = zeros(2,p.nsims)

	N = Normal(0,1)  # normal dist N(0,σ)

	# random setup
	Random.seed!(p.rseed)
	wshocks = quantile.(N, rand(p.nsims,p.nT)) .* p.sigma
	dshocks = rand(p.nsims,p.nT) # discrete choice shocks
	rshocks = rand(p.nsims,p.nT) # returning to work shocs
	w0shocks = rand(p.nsims)



	for it in 1:p.nT
		if it == 1
			# draw from initial distributions
			s.w0[:, it] = p.initw0 .+ w0shocks *(p.initw1 - p.initw0)
			s.worker[:, it] .= 1   # all work in frist period
		else

			# working state
			s.worker[ working, it] .= 1 # chose working, so stay worker
			s.worker[.!(working), it] .= 2 # chose retire
			s.ret_age[(s.worker[:,it-1] .== 1) .& (s.worker[:,it] .== 2)] .= it

			# income
			s.inc[.!(working), it] .= p.pension
			s.inc[ working, it] .= income(it,p,wshocks[working,it])


			# end of period wealth in period it
			s.w0[ :     , it] .= s.w1[:      , it-1]*p.R
			s.w0[working, it] .= s.w1[working, it-1]*p.R .+ s.inc[ working, it]
		end
		# consumption
		s.cons[ working,it]    = gety(interp(m.c[1,it].env,s.w0[ working, it]))

		s.cons[.!(working),it] = gety(interp(m.c[2,it].env,s.w0[.!(working), it]))
		# end of period wealth
		s.w1[:,it] = s.w0[:,it] - s.cons[:,it]
		clamp!(s.w1,p.a_low, p.a_high)

		if any(s.cons[ :,it] .< 0)
			println("neg cons")
			println(s.cons[s.cons[ :,it] .< 0, it])
			println(s.w0[.!(working), it])
			error()
		end

		# CCP to remain worker
		vmat[1,:] = vfun(1,it, s.cons[:,it] ,s.w0[:,it], m.v[1,it],p )
		vmat[2,:] = vfun(2,it, s.cons[:,it] ,s.w0[:,it], m.v[2,it],p )
		s.prob_work[:,it] = (s.worker[:,it] .== 1) .* ccp(vmat,p)  # prob to remain worker

		# discrete choice of workers
		working[:] = s.prob_work[:,it] .> dshocks[:,it]

		# random transition back into work state for retirees
		working[.!(working)] = p.delta .> rshocks[.!(working),it]
	end
	s
end


function sim(m::GModel,p::Param)

	s = Simulation(p)

	# random
	Random.seed!(p.rseed)
	wshocks = rand(p.nsims,p.nT)
	dshocks = rand(p.nsims,p.nT)
	w0shocks = rand(p.nsims)
	s.ystate[:,1] .= rand(1:p.ny,p.nsims)  # randomly allocate to an income state

	iy = 0 # y-state index
	working   = trues(p.nsims)  # indicator of currentlcy working or not
	vmat      = zeros(2,p.nsims)
	# s.shocks[:] = rand(p.nT * p.nsims)
	GG = cumsum(m.ywgt,dims = 2)
	for it in 2:p.nT
		for i in 1:p.nsims
			s.ystate[i,it] = searchsortedfirst(GG[s.ystate[i,it-1],:], wshocks[i,it])
		end
	end

	for it in 1:p.nT
		if it == 1
			# draw from initial distributions
			s.w0[:, it] = p.initw0 .+ w0shocks * (p.initw1 - p.initw0)
			s.worker[:, it] .= 1   # all work in frist period
			s.inc[ working, it] .= income(it,p,m.yvec[s.ystate[working,it]])
		else

			# working state
			s.worker[ working, it] .= 1 # chose working, so stay worker
			s.worker[.!(working), it] .= 2 # chose retire
			s.ret_age[(s.worker[:,it-1] .== 1) .& (s.worker[:,it] .== 2)] .= it

			# income
			s.inc[.!(working), it] .= 0.0
			s.inc[ working, it] .= income(it,p,m.yvec[s.ystate[working,it]])


			# start of period wealth in period it
			s.w0[ :     , it] .= s.w1[:      , it-1]*p.R
			s.w0[working, it] .= s.w1[working, it-1]*p.R .+ s.inc[ working, it]
		end
		# consumption and value
		# can be made much faster by grouping individuals by `iy` state
		for i in 1:p.nsims
			iy = s.ystate[i,it]
			if working[i]
				s.cons[i,it]    = gety(interp(m.c[1,iy,it].env,[s.w0[ i, it]]))[1]
			else
				s.cons[i,it]    = gety(interp(m.c[2,iy,it].env,[s.w0[ i, it]]))[1]
			end
			vmat[1,i] = vfun(1,it, [s.cons[i,it] ],[s.w0[i,it]], m.v[1,iy,it],p )[1]
			vmat[2,i] = vfun(2,it, [s.cons[i,it] ],[s.w0[i,it]], m.v[2,iy,it],p )[1]
		end
		# end of period wealth
		s.w1[:,it] = s.w0[:,it] - s.cons[:,it]

		# CCP to remain worker
		s.prob_work[:,it] = (s.worker[:,it] .== 1) .* ccp(vmat,p)  # prob to remain worker

		# discrete choice
		working[:] = s.prob_work[:,it] .> dshocks[:,it]
	end
	s
end

function sim(m::BModel,p::Param)

	s = BSimulation(p)

	# s.bkflag: 1 = flag off, 2 = flag on

	# random
	Random.seed!(p.rseed)
	wshocks = rand(p.nsims,p.nT)
	dshocks = rand(p.nsims,p.nT)
	w0shocks = rand(p.nsims)
	s.ystate[:,1] .= rand(1:p.ny,p.nsims)  # randomly allocate to an income state

	iy = 0 # y-state index

	# everybody starts as not
	filing    = falses(p.nsims)  # indicator of whether is not filing in current period
	notfiling = trues(p.nsims)
	vmat      = zeros(2,p.nsims)

	# s.shocks[:] = rand(p.nT * p.nsims)
	GG = cumsum(m.ywgt,dims = 2)
	for it in 2:p.nT
		for i in 1:p.nsims
			s.ystate[i,it] = searchsortedfirst(GG[s.ystate[i,it-1],:], wshocks[i,it])
		end
	end

	for it in 1:p.nT
		# println(it)
		notfiling[:] .= .!(filing)
		if it == 1
			# draw from initial distributions
			s.w0[:, it] = p.initw0 .+ w0shocks * (p.initw1 - p.initw0)
			s.ainit[:,it] = s.w0[:, it]

			# everybody starts with bkflag off
			s.inc[: , it] .= income(it,p,m.yvec[s.ystate[:,it]])  # all get their income
		else

			# set bk flag at start of period depending on last periods choices:
			# if not filing, flag does not change
			s.bkflag[notfiling , it] = s.bkflag[notfiling , it-1]
			# else, flag goes from false to true
			s.bkflag[filing , it] .= true  # have chosen file in previous period, so flag is on
			s.file_age[filing] .= it  # age at which filed

			# income
			s.inc[: , it] .= income(it,p,m.yvec[s.ystate[:,it]])  # all get their income

			# s.inc[.!(notfiling), it] .= 0.0
			# s.inc[ notfiling, it] .= income(it,p,m.yvec[s.ystate[notfiling,it]])


			# start of period wealth in period it
			# s.w0[ :     , it] .= s.w1[:      , it-1]*p.R
			s.ainit[:,it] = s.w1[:, it-1]*p.R
			s.w0[notfiling, it] .= s.ainit[notfiling,it] .+ s.inc[notfiling, it]
			s.w0[filing    , it] .= 0.0 .+ s.inc[filing, it]
		end

		# end of recording current state
		fill!(filing, false)

		# figure out discrete choice in non bk state first
		inoflag = findall(.!(s.bkflag[:,it]))
		  iflag = findall(   s.bkflag[:,it])
		for i in inoflag
			iy = s.ystate[i,it]
			vmat[1,i] = gety(interp(m.v[1,iy,it].env,[s.w0[ i, it]]))[1]
			vmat[2,i] = s.ainit[i,it] < 0 ? m.v[2,iy,it].env[1].y : -Inf   # file only with neg assets
			s.vnofile[i,it] = vmat[1,i]
			s.vfile[i,it] = vmat[2,i]
		end


		s.prob_nofile[inoflag,it] = ccp(vmat[:,inoflag],p)  # prob to not file
		s.prob_nofile[iflag  ,it] .= NaN   # no discrete choice with a flag on
		# discrete choice
		filing[inoflag] = s.prob_nofile[inoflag,it] .< dshocks[inoflag,it]
		filing[iflag]  .= false

		# consumption depends on whehter flag and if not, which dchoice
		inoflag_nofile = findall( .!(s.bkflag[:,it]) .& .!(filing) )
		inoflag_file   = findall( .!(s.bkflag[:,it]) .& filing )
		for i in inoflag_nofile
			s.cons[i,it] = gety(interp(m.c[1,s.ystate[i,it],it].env,[s.w0[ i, it]]))[1]
			# end of period wealth
			s.w1[i,it] = s.w0[i,it] - s.cons[i,it]
		end
		for i in inoflag_file
			s.cons[i,it]    = m.c[2,s.ystate[i,it],it].env[1].y
			# end of period wealth
			s.w1[i,it] = 0.0
		end
		for i in iflag
			s.cons[i,it] = gety(interp(m.cbk[s.ystate[i,it],it].env,[s.w0[ i, it]]))[1]
			s.w1[i,it] = s.w0[i,it] - s.cons[i,it]

		end

	end
	s
end

function rbs()
	m,p = DCEGM.runbk()
	sim(m,p)
end

function rs(;par = Dict())
	m,p = DCEGM.runf(par = par)
	sim(m,p)
end
function rsp(;par = Dict())
	plot_s(rs(par=par))
end

function rsg(;par = Dict())
	m,p = DCEGM.rung(par = par)
	sim(m,p)
end
function rsgp(;par = Dict())
	plot_s(rsg(par=par))
end

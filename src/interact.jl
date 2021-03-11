
function iminimal()
	p = Param()
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	sigmas = 0.0:0.05:1.0

	mp = @manipulate for γ in slider(gammas, label = "γ", value =p.gamma ),
						 β in slider(betas, label = "β", value =p.beta) ,
						 R in slider(Rs, label = "R", value =p.R) ,
						 σ in slider(sigmas, label = "σ", value =p.sigma)

		 p =Param(par = Dict(:beta => β, :sigma => σ,  :R => R))
		 minimal_EGM(p)

	end
	@layout! mp vbox(
		hbox(β, R, σ),
		observe(_))
	w = Window()
	body!(w,mp)
end


function iminimalb()
	p = Param()
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	sigmas = 0.0:0.05:1.0
	nus = 0.0:0.01:10
	bbars = 0.0:0.1:10

	mp = @manipulate for γ in slider(gammas, label = "γ", value =p.gamma ),
						 β in slider(betas, label = "β", value =p.beta) ,
						 R in slider(Rs, label = "R", value =p.R) ,
						 σ in slider(sigmas, label = "σ", value =p.sigma),
						 bbar in slider(bbars, label = "bbar", value =p.bbar),
						 nu in slider(nus, label = "ν", value =5.0)

		 p =Param(par = Dict(:beta => β, :sigma => σ,  :R => R , :bbar => bbar, :ν => nu))
		 minimal_EGM_bequest(p)

	end
	@layout! mp vbox(
		hbox(β, R, σ),
		hbox(bbar, nu),
		observe(_))
	w = Window()
	body!(w,mp)
end


# DCEGM.interact(DCEGM.rsgp) for sim
# DCEGM.interact(DCEGM.runfp) for fedors model
function ifedor()
	p = Param()
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	alphas = 0.0:0.05:0.5
	sigmas = 0.0:0.05:0.35
	lambdas = 0.0000002:0.05:1
	rhos = 0.1:0.05:1
	deltas = 0.0:0.05:0.5
	pensions = 0.0:0.1:1.0
	nus = 0.0:0.1:10
	bbars = 0.0:0.1:10

	mp = @manipulate for dosim = Dict("sim" => true, "sol" => false),
						 id = Dict("id=$id" => id for id in 1:2) ,
						 nsims = spinbox(label="nsims"; value=p.nsims) |> onchange,
						 # γ in slider(gammas, label = "γ", value =p.gamma ),
						 β in slider(betas, label = "β", value =p.beta) |> onchange ,
						 R in slider(Rs, label = "R", value =p.R) |> onchange,
						 α in slider(alphas, label = "α", value =p.alpha) |> onchange ,
						 σ in slider(sigmas, label = "σ", value =p.sigma) |> onchange ,
						 λ in slider(lambdas, label = "λ", value =p.lambda) |> onchange,
						 δ in slider(deltas, label= "δ",value = p.delta) |> onchange,
						 pens in slider(pensions, label= "pension",value = p.pension) |> onchange,
						 bbar in slider(bbars, label = "bbar", value =p.bbar) |> onchange,
						 nu in slider(nus, label = "ν", value =p.ν) |> onchange

		 m,p = runf(par = Dict(:nsims => nsims , :bbar => bbar, :ν => nu, :beta => β, :alpha => α, :sigma => σ, :lambda => λ, :R => R, :delta => δ, :pension => pens))
		 if dosim
			 s = sim(m,p)
			 plot_s(s)
		 else
			 plot(m,p,id = id)
		 end
	end
	@layout! mp vbox(
		hbox(dosim, id, nsims),
		hbox(β, R, σ),
		hbox(α, λ, nu),
		hbox(δ, pens, bbar),
		observe(_))
	w = Window()
	body!(w,mp)
end

function igmodel()
	p = Param()
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	alphas = 0.0:0.05:0.5
	sigmas = 0.0:0.05:0.35
	lambdas = 0.0000002:0.05:1
	rhos = 0.1:0.05:1
	deltas = 0.0:0.05:0.5
	pensions = 0.0:0.1:1.0
	nus = 0.0:0.1:10
	bbars = 0.0:0.1:10

	try
		mp = @manipulate for dosim = Dict("sim" => true, "sol" => false),
							 id = Dict("id=$id" => id for id in 1:2),
							 nsims = spinbox(label="nsims"; value=20) |> onchange,
							 # γ in slider(gammas, label = "γ", value =p.gamma ),
							 β in slider(betas, label = "β", value =p.beta) |> onchange ,
							 R in slider(Rs, label = "R", value =p.R)  |> onchange,
							 α in slider(alphas, label = "α", value =p.alpha) |> onchange ,
							 σ in slider(sigmas, label = "σ", value =p.sigma) |> onchange ,
							 λ in slider(lambdas, label = "λ", value =p.lambda) |> onchange,
							 ρ in slider(0:0.05:1, label = "ρ", value =p.ρ) |> onchange,
							 δ in slider(deltas, label= "δ",value = p.delta) |> onchange,
							 pens in slider(pensions, label= "pension",value = p.pension) |> onchange,
							 bbar in slider(bbars, label = "bbar", value =p.bbar) |> onchange,
							 nu in slider(nus, label = "ν", value =p.ν) |> onchange

			 m,p = rung(par = Dict(:nsims => nsims , :bbar => bbar, :ν => nu, :ρ => ρ, :beta => β, :alpha => α, :sigma => σ, :lambda => λ, :R => R, :delta => δ, :pension => pens))
			 if dosim
				 s = sim(m,p)
				 plot_s(s)
			 else
				 plot(m,p,id = id,iy = 1)
			 end
		end
		@layout! mp vbox(
			hbox(dosim, id, nsims),
			hbox(β, R, σ),
			hbox(α, λ, nu, ρ),
			hbox(δ, pens, bbar),
			observe(_))
		w = Window()
		body!(w,mp)
	catch e
		return e
	end

end

function ibkmodel(;it::Bool=false)
	p = Param(par = Dict(:nT => 20, :retage => 10) )
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	alphas = [0.0,0.1,0.35,0.5,1.0]
	alphaT = -1000:1:0.0
	lambdas = 0.2:0.1:1.0
	rhos = 0.1:0.05:1
	times = 1:p.nT

	if it
		mp = @manipulate for id = Dict("id=$id" => id for id in 1:2),
							 ti in slider(times, label = "period", value = 1),
			                 γ in slider(gammas, label = "γ", value =p.gamma ),
							 β in slider(betas, label = "β", value =p.beta) ,
							 a in slider(alphas, label = "α", value = p.alpha) ,
							 aT in slider(alphaT, label = "αT", value = p.alphaT) ,
							 alowT in slider(-5:0.5:0.0, label = "alow") ,
							 λ in slider(lambdas, label = "λ", value =p.lambda)#,
							 # ρ in slider(rhos, label = "ρ", value =p.ρ)

			m,p = runbk(par = Dict(:nT => 50,:a_low => -5.0,:a_lowT => alowT,:na =>501,:beta => β, :alphaT => aT, :alpha => a, :gamma => γ , :lambda => λ))
			plot(m,p,it = ti, id = id,ylims = (-5,13),size = (700,400))
		end
		w = Window()
		body!(w,mp)

	else
		mp = @manipulate for iy = OrderedDict("iy=$iy" => iy for iy in 1:p.ny) ,
			                 id = Dict("id=$id" => id for id in 1:2),
			                 γ in slider(gammas, label = "γ", value =p.gamma ),
							 β in slider(betas, label = "β", value =p.beta) ,
							 a in slider(alphas, label = "α", value = p.alpha) ,
							 aT in slider(alphaT, label = "αT", value = p.alphaT) ,
							 alow in slider(-5:0.5:0.0, label = "alow") ,
							 λ in slider(lambdas, label = "λ", value =p.lambda)#,
							 # ρ in slider(rhos, label = "ρ", value =p.ρ)

			m,p = runbk(par = Dict(:nT => 50,:a_low => -5.0,:a_lowT => alow,:na =>501,:beta => β, :alphaT => aT, :alpha => a, :gamma => γ , :lambda => λ))
			plot(m,p,iy = iy, id = id,ylims = (-5,13),size = (700,400))
		end
		w = Window()
		body!(w,mp)

	end


end


function ibksim(;pars = Dict(:nT => 50),it::Bool=false)
	p = Param(par = pars )
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	alphas = [0.0,0.1,0.35,0.5,1.0]
	alphaT = -10:0.5:0.0
	lambdas = 0.02:0.1:1.0
	rhos = 0.1:0.05:1
	times = 1:p.nT

	if it
		mp = @manipulate for dosim = Dict("sim" => true, "sol" => false),
			id = Dict("id=$id" => id for id in 1:2),
							 ti in slider(times, label = "period", value = 1),
			                 γ in slider(gammas, label = "γ", value =p.gamma ),
							 β in slider(betas, label = "β", value =p.beta) ,
							 a in slider(alphas, label = "α", value = p.alpha) ,
							 aT in slider(alphaT, label = "αT", value = p.alphaT) ,
							 alowT in slider(-5:0.5:0.0, label = "alowT") ,
							 λ in slider(lambdas, label = "λ", value =p.lambda)#,
							 # ρ in slider(rhos, label = "ρ", value =p.ρ)
			pp = merge(pars,Dict(:a_low => -5.0,:a_lowT => alowT,:na =>501,:beta => β, :alphaT => aT, :alpha => a, :gamma => γ , :lambda => λ))
			m,p = runbk(par = pp)
			if dosim
				s = sim(m,p)
				plot_s(s)
			else
				plot(m,p,it = ti, id = id,ylims = (-5,13),size = (700,400))
			end

		end
		w = Window()
		body!(w,mp)

	else
		mp = @manipulate for dosim = Dict("sim" => true, "sol" => false),
			iy = OrderedDict("iy=$iy" => iy for iy in 1:p.ny) ,
			                 id = Dict("id=$id" => id for id in 1:2),
			                 γ in slider(gammas, label = "γ", value =p.gamma ),
							 β in slider(betas, label = "β", value =p.beta) ,
							 a in slider(alphas, label = "α", value = p.alpha) ,
							 aT in slider(alphaT, label = "αT", value = p.alphaT) ,
							 alowT in slider(-5:0.5:0.0, label = "alowT", value = p.a_lowT) ,
							 λ in slider(lambdas, label = "λ", value =p.lambda)#,
							 # ρ in slider(rhos, label = "ρ", value =p.ρ)

			pp = merge(pars,Dict(:a_low => -5.0,:a_lowT => alowT,:na =>101,:beta => β, :alphaT => aT, :alpha => a, :gamma => γ , :lambda => λ))
 			m,p = runbk(par = pp)
			if dosim
 				s = sim(m,p)
 				plot_s(s)
 			else
 				plot(m,p,iy = iy, id = id,ylims = (-5,13),size = (700,400))
 			end
		end
		w = Window()
		body!(w,mp)
	end
end



# DCEGM.interact(DCEGM.rsgp) for sim
# DCEGM.interact(DCEGM.runfp) for fedors model
function interact(fun::Function)
	p = Param()
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	alphas = [0.0,0.35,0.5]
	sigmas = 0.0:0.05:0.35
	lambdas = 0.0000002:0.05:1
	rhos = 0.1:0.05:1

	mp = @manipulate for γ in slider(gammas, label = "γ", value =p.gamma ),
						 β in slider(betas, label = "β", value =p.beta) ,
						 R in slider(Rs, label = "R", value =p.R) ,
						 α in slider(alphas, label = "α", value =p.alpha) ,
						 σ in slider(sigmas, label = "σ", value =p.sigma) ,
						 λ in slider(lambdas, label = "λ", value =p.lambda),
						 ρ in slider(rhos, label = "ρ", value =p.ρ)

		fun(par = Dict(:ρ => ρ, :gamma => γ, :beta => β, :alpha => α, :sigma => σ, :lambda => λ, :R => R))
	end
end

function igmodel()
	p = Param()
	gammas = 1.0:0.1:3.0
	betas  = 0.5:0.05:1.0
	Rs  = 1.0:0.05:1.5
	alphas = [0.0,0.35,0.5]
	sigmas = 0.0:0.05:0.35
	lambdas = 0.0000002:0.05:1
	rhos = 0.1:0.05:1

	mp = @manipulate for iy = OrderedDict("iy=$iy" => iy for iy in 1:p.ny) ,
		                 id = Dict("id=$id" => id for id in 1:2),
		                 γ in slider(gammas, label = "γ", value =p.gamma ),
						 β in slider(betas, label = "β", value =p.beta) ,
						 σ in slider(sigmas, label = "σ", value =p.sigma) ,
						 λ in slider(lambdas, label = "λ", value =p.lambda),
						 ρ in slider(rhos, label = "ρ", value =p.ρ)

		m,p = rung(par = Dict(:ρ => ρ, :gamma => γ, :beta => β, :sigma => σ, :lambda => λ))
		plot(m,p,iy = iy, id = id,ylims = (-15,13),size = (700,400))
	end
end

function ibkmodel(;it::Bool=false)
	p = Param()
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
							 λ in slider(lambdas, label = "λ", value =p.lambda)#,
							 # ρ in slider(rhos, label = "ρ", value =p.ρ)

			m,p = runbk(par = Dict(:a_low => -5.0,:a_lowT => -5.0,:na =>101,:beta => β, :alphaT => aT, :alpha => a, :gamma => γ , :lambda => λ))
			plot(m,p,it = ti, id = id,ylims = (-5,13),size = (700,400))
		end

	else
		mp = @manipulate for iy = OrderedDict("iy=$iy" => iy for iy in 1:p.ny) ,
			                 id = Dict("id=$id" => id for id in 1:2),
			                 γ in slider(gammas, label = "γ", value =p.gamma ),
							 β in slider(betas, label = "β", value =p.beta) ,
							 a in slider(alphas, label = "α", value = p.alpha) ,
							 aT in slider(alphaT, label = "αT", value = p.alphaT) ,
							 λ in slider(lambdas, label = "λ", value =p.lambda)#,
							 # ρ in slider(rhos, label = "ρ", value =p.ρ)

			m,p = runbk(par = Dict(:a_low => -5.0,:a_lowT => -5.0,:na =>101,:beta => β, :alphaT => aT, :alpha => a, :gamma => γ , :lambda => λ))
			plot(m,p,iy = iy, id = id,ylims = (-5,13),size = (700,400))
		end

	end


end

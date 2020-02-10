


@testset "F(edor) Test" begin

	# set up git repo. that's a submodule at matlab/
	tdir = @__DIR__
	mldir = joinpath(tdir,"..","matlab")
	cd(mldir)
	br = chomp(read(`git rev-parse --abbrev-ref HEAD`,String))
	if br !="bm"
		run(`git checkout bm`)
	end

	# run fedor's matlab code
	# if haskey(ENV,"TRAVIS")
	# 	println("no matlab license on travis.")
	# 	println("will use saved results instead")
	# else
	# 	run(`/Applications/MATLAB_R2019b.app/bin/matlab -batch "bench"`)
	# end

	cd(tdir)

	# set up julia model with identical parameter settings:
	# m5=model_retirement;
	# m5.ngridm=500;
	# m5.df=1/(1+m5.r); %flat consumption hopefully
	# m5.sigma=0.35;
	# m5.lambda=0.000002;

	# run benchmark
	# 1. warm up julia JIT on a small version
	pd = Dict(:na => 500,
			   :beta => 0.95,
			   :sigma => 0.35,
			   :R => 1.05,
			   :lambda => 0.000002)

	# DCEGM.runf(par = pd)
	pd[:na] = 500
	m,p = DCEGM.runf(par = pd)

	# read matlab results and test against each julia result set
	for i in 1:p.nD
		for it in p.nT:-1:1
			mpo = readdlm(joinpath(mldir,"output","policy_$(i)_$(it).csv"), ',')
			mvf = readdlm(joinpath(mldir,"output","value_$(i)_$(it).csv"), ',')

			@testset "period $it choice $i" begin
				@test length( getx(m.c[i,it].env) ) == length( mpo[1,:] )
				@test length( getx(m.v[i,it].env) ) == length( mvf[1,:] )
				fedsort = sortperm(mpo[1,:]) # sometimes fedors policy is not sorted!
				@test maximum( abs.( getx(m.c[i,it].env) .- mpo[1,fedsort] ) )  < 1e-4
				@test maximum( abs.( gety(m.c[i,it].env) .- mpo[2,fedsort] ) )  < 1e-4

				if it==p.nT
					@test all(isapprox.( getx(m.v[i,it].env)[1] , mvf[1,1] , atol=1e-4) )
					@test all(isapprox.( gety(m.v[i,it].env)[1] , mvf[2,1] , atol=1e-4) )
				else
					@test maximum( abs.( getx(m.v[i,it].env) .- mvf[1,:] ) ) < 1e-4
					@test maximum( abs.( gety(m.v[i,it].env) .- mvf[2,:] ) ) < 5e-4  # numerical imprecision from ACII export
				end
			end

		end
	end
end

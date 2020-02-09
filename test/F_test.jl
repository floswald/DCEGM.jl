


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
	if haskey(ENV,"TRAVIS")
		println("no matlab license on travis.")
		println("will use saved results instead")
	else 
		run(`/Applications/MATLAB_R2019b.app/bin/matlab -batch "bench"`)
	end

	cd(tdir)

	# set up julia model with identical parameter settings:
	# m5=model_retirement;
	# m5.ngridm=500;
	# m5.df=1/(1+m5.r); %flat consumption hopefully
	# m5.sigma=0.35;
	# m5.lambda=0.000002;

	# run benchmark
	# 1. warm up julia JIT on a small version
	pd = Dict(:na => 15, 
			   :beta => 0.95, 
			   :sigma => 0.35,
			   :R => 1.05,
			   :lambda => 0.000002)

	DCEGM.runf(par = pd)

	# 2. take actual timing
	pd[:na] = 500
	m,p = @time DCEGM.runf(par = pd)

	# read matlab results and test against each julia result set
	for i in 1:p.nD
		for it in 1:p.nT
			mpo = readdlm(joinpath(mldir,"output","policy_$(i)_$(it).csv"), ',')
			mvf = readdlm(joinpath(mldir,"output","value_$(i)_$(it).csv"), ',')

			@test getx(m.c[i,it]) .≈ mpo[1,:]
			@test gety(m.c[i,it]) .≈ mpo[2,:]

			@test getx(m.v[i,it]) .≈ mvf[1,:]
			@test gety(m.v[i,it]) .≈ mvf[2,:]

		end
	end
end
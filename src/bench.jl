

function bm()
	# set up git repo. that's a submodule at matlab/
	tdir = @__DIR__
	mldir = joinpath(tdir,"..","matlab")
	cd(mldir)
	br = chomp(read(`git rev-parse --abbrev-ref HEAD`,String))
	if br !="bm"
		run(`git checkout bm`)
	end

	# run fedor's matlab code
	run(`/Applications/MATLAB_R2019b.app/bin/matlab -batch "bench"`)

	cd(joinpath(tdir,".."))
	# set up param dict
    pd = Dict(:na => 500,
			   :beta => 0.95,
			   :sigma => 0.35,
			   :R => 1.05,
			   :lambda => 0.000002)

	# accuracy of solution has been tested in F_test unit test.
	println()
    DCEGM.runf(par = pd);  # run once to compile
	println("julia timing:")
	@benchmark DCEGM.runf(par = $pd)
	return nothing

end

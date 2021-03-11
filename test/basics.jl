
@testset "deltatest" begin
	p=DCEGM.Param(par=Dict(:delta => 0.9))
	@test p.delta == 0.9
end

@testset "basics" begin
	p = DCEGM.Param()
	@testset "run minimal egm" begin
		m = DCEGM.minimal_EGM(p)
		@test true
	end

	@testset "ccp" begin
		x = rand(5,10)
		c = DCEGM.ccp(x,p)
		@test length(c) == 10
		@test minimum(c) >= 0.0
		@test maximum(c) <= 1.0
	end
	@testset "logsum" begin
		x = rand(5,10)
		c = DCEGM.logsum(x,p)
		@test length(c) == 10
		@test maximum(abs,DCEGM.logsum(zeros(2,2),p) .- p.lambda * [log(2) log(2)]) .< 1e-6
		p.lambda = 1
		@test maximum(abs,DCEGM.logsum(zeros(2,2),p) .- p.lambda * [log(2) log(2)]) .< 1e-6
	end
end

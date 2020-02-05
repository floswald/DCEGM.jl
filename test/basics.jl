



@testset "components" begin
	p = DCEGM.Param()
	@testset "run minimal egm" begin
		m = DCEGM.minimal_EGM()
		@test m[1][end][1] == 0.0
		@test m[1][end][end] == p.a_high
		@test m[2][end][1] == 0.0
		@test m[2][end][end] == p.a_high
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
	# @testset "borrowing constraints" begin
	# m,p = DCEGM.model()
	#     @test DCEGM.nextbound(m.yvec[1],p.nT-2,p) == -DCEGM.income(p.nT-1,p,m.yvec[1])*p.R^(-1)
	#     @test DCEGM.nextbound(m.yvec[1],p.nT-3,p) == -DCEGM.income(p.nT-2,p,m.yvec[1])*p.R^(-1) - DCEGM.income(p.nT-1,p,m.yvec[1])*p.R^(-2)
	# end
end

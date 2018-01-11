using DCEGM
using Base.Test

# write your own tests here
@testset "run minimal egm" begin
	p = DCEGM.Param()
	m = DCEGM.minimal_EGM()
	@test m[1][end][1] == p.a_low
	@test m[1][end][end] == p.a_high
	@test m[2][end][1] == p.a_low
	@test m[2][end][end] == p.a_high
end


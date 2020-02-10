

@testset "Envelope" begin
    @testset "Constructors" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        L1 = MLine(x1,x1)

        en = Envelope(L1)
        @test isa(en,Envelope)
        @test length(en.env)==n
        @test length(gets(en))==0
        DCEGM.removed!(en,L1)
        @test length(getr(en))==0
        @test en.vbound == 0.0
    end
end

@testset "Testing upper_env" begin
    @testset "upper_env test 1" begin

        x1,x2,a = DCEGM.tplot1()
        e = DCEGM.upper_env(a)
        @test issorted(getx(e.env))

        # test actual values in the Envelope
        @test all( gety(filter(x -> x.x < 5.0 , e.env.v)) .== 5.0 )
        @test all( gety(filter(x -> x.x > 5.0 , e.env.v)) .== getx(filter(x -> x.x > 5.0 , e.env.v)) )

    end
    # @testset "upper_env test 2" begin
    #     n = 15
    #     x1 = collect(range(0,stop = 10,length = n))
    #     L1 = MLine(x1,x1)
    #     L2 = MLine(x1,ones(n)*5)
    #     e = Envelope([L1,L2])
    #     upper_env!(e)
    #     @test issorted(getx(e.env))
    #     @test all( gety(filter(x -> x.x < 5.0 , e.env.v)) .== 5.0 )
    #     @test all( gety(filter(x -> x.x > 5.0 , e.env.v)) .== getx(filter(x -> x.x > 5.0 , e.env.v)) )
    #
    #     upper_env!(e)
    #     @test issorted(getx(e.env))
    #     @test gets(e)[1].x ≈ 5.0
    #     @test gets(e)[1].y ≈ 5.0
    # end
    # @testset "upper_env test 3" begin
    #     n = 15
    #     x1 = collect(range(0,stop = 10,length = n))
    #     L1 = MLine(x1,x1)
    #     L2 = MLine(x1,5.0 .+ 0.3*x1)
    #     e = Envelope([L1,L2])
    #     upper_env!(e)
    #     @test issorted(getx(e.env))
    #     @test all( gety(filter(x -> x.x < 7.0 , e.env.v)) .== gety(filter(x -> x.x < 7.0 , L2.v)) )
    #
    #     upper_env!(e)
    #     @test issorted(getx(e.env))
    #     @test gets(e)[1].x ≈ 5.0/0.7
    #     @test gets(e)[1].y ≈ 5.0/0.7
    # end
    # @testset "upper_env test 4" begin
    #     n = 15
    #     x1 = collect(range(0,stop = 10,length = n))
    #     x2 = collect(range(-1,stop = 9,length = n))
    #     L1 = MLine(x1,x1)
    #     L2 = MLine(x2,5.0 .+ 0.3*x2)
    #     e = Envelope([L1,L2])
    #     upper_env!(e)
    #     @test issorted(getx(e.env))
    #     # @test isapprox(getx(e.env),sort(unique(vcat(x1,x2))), atol = 1e-15)
    #     @test gets(e)[1].x ≈ 5.0/0.7
    #     @test gets(e)[1].y ≈ 5.0/0.7
    # end
#     # poor test since we only do linear interpolation...
#     @testset "log(x) vs 2(x-1)" begin
#         x = collect(0.0:0.1:1.5)
#         L1 = MLine(x,log.(x))
#         L2 = MLine(x, 2 .* (x .- 1.0) )
#         e = Envelope([L1,L2])
#         upper_env!(e)
#         @test issorted(getx(e.env))
#         # @test gets(e)[1].x ≈ DCEGM.fzero(x-> 2*(x-1) - log(x),0.01,0.5)
#         # @test gets(e)[1].y ≈ log( gets(e)[1].x )
#         @test gets(e)[2].x == 1.0
#         @test gets(e)[2].y == 0.0
#     end
    @testset "y=x, y=2x, y=0.1x + 1.98" begin
        x1 = collect(-0.9:0.3:2.7)
        L1 = MLine(x1, x1)
        x2 = collect(0.0:0.1:1)
        L2 = MLine(x2, 2 .* x2)
        x3 = collect(1.0:0.45:2.9)
        L3 = MLine(x3, (0.1 .* x3) .+ 1.9)
        a = [L1,L2,L3]
        e = DCEGM.upper_env(a)
        @test issorted(getx(e.env))
        @test all( gety(filter(x -> x.x < 0.0 , e.env.v)) .== gety(filter(x -> x.x < 0.0 , L1.v)) )
        @test all( gety(filter(x -> (x.x > 0.0) & (x.x < 1.0) , e.env.v)) .== gety(filter(x -> (x.x > 0.0) & (x.x < 1.0)  , L2.v)) )

        @test gets(e)[1] ≈ Point(1.9/0.9,1.9/0.9)

    end
    @testset "upper_env test: decreasing " begin
        L,e,x1,x2,pl = DCEGM.test_upper_env_dec()
        # n = 15
        # x1 = collect(range(0,stop = 10,length = n))
        # x2 = collect(range(-1,stop = 9,length = n))
        # x = vcat(x1,x2)
        # y = vcat(x1[end:-1:1],ones(n)*5)
        # L = MLine(x,y)
        # e = splitLine(L)
        # upper_env!(e)
        @test issorted(getx(e.env))
    end
end
@testset "splitLine" begin
    @testset "test 1: simple" begin
        x = [1,2,3,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.2,1.8,2.1]
        L1 = MLine(x,y)
        e = splitLine(L1)

        @test isa(e,Vector{MLine{Float64}})
        @test length(e) == 3
    end
    @testset "test 2: less simple" begin
        x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
        L1 = MLine(x,y)
        a = splitLine(L1)
        @test isa(a,Array{MLine{Float64}})
        e = DCEGM.upper_env(a)
        @test isa(e,Envelope)
        @test length(gets(e))==1
        @test length(a) == 3
        @test eltype(a)== MLine{Float64}

        @test extrema(getx(a[1].v)) == (1.0,3.0)
        @test extrema(getx(a[2].v)) == (1.5,3.0)
        @test extrema(getx(a[3].v)) == (1.5,2.9)
    end
end

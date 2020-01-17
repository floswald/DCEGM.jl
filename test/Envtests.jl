

@testset "Envelope" begin
    @testset "Constructors" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        L1 = MLine(x1,x1)

        en = Envelope(L1)
        @test isa(en,Envelope)
        @test length(en.L)==0
        @test length(en.env)==n
        @test length(gets(en))==0
        @test length(getr(en))==0
        DCEGM.removed!(en)
        @test length(getr(en))==0
        @test en.env_set
        @test en.vbound == 0.0

        en = Envelope([L1,L1])
        @test isa(en,Envelope)
        @test length(en.L)==2
        @test length(en.env)==1
        @test length(gets(en))==0
        @test length(getr(en))==0
        @test !en.env_set

        en = Envelope(1)
        @test isa(en,Envelope{Int})
        @test !en.env_set
        @test length(en.env)==1
    end
end

@testset "Testing Envelopes over 2 MLines" begin
    @testset "upper_env test 1" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        x2 = collect(range(-1,stop = 9,length = n))
        L1 = MLine(x1,x1)
        L2 = MLine(x2,ones(n)*5)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test issorted(getx(e.env))
        @test isapprox(getx(e.env),sort(unique(vcat(x1,x2))), atol = 1e-15)

        # test actual values in the Envelope
        @test all( gety(filter(x -> x.x < 5.0 , e.env.v)) .== 5.0 )
        @test all( gety(filter(x -> x.x > 5.0 , e.env.v)) .== getx(filter(x -> x.x > 5.0 , e.env.v)) )

        upper_env!(e, do_intersect = true)
        @test issorted(getx(e.env))
        @test isapprox(getx(e.env),sort(unique(vcat(x1,x2))), atol = 1e-15)
        @test gets(e)[1].x == 5.0
        @test gets(e)[1].y == 5.0
    end
    @testset "upper_env test 2" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        L1 = MLine(x1,x1)
        L2 = MLine(x1,ones(n)*5)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test getx(e.env) == x1
        @test issorted(getx(e.env))
        @test all( gety(filter(x -> x.x < 5.0 , e.env.v)) .== 5.0 )
        @test all( gety(filter(x -> x.x > 5.0 , e.env.v)) .== getx(filter(x -> x.x > 5.0 , e.env.v)) )

        upper_env!(e, do_intersect = true)
        @test getx(e.env) == x1
        @test issorted(getx(e.env))
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
    end
    @testset "upper_env test 3" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        L1 = MLine(x1,x1)
        L2 = MLine(x1,5.0 .+ 0.3*x1)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test isapprox(getx(e.env),x1)
        @test issorted(getx(e.env))
        @test all( gety(filter(x -> x.x < 7.0 , e.env.v)) .== gety(filter(x -> x.x < 7.0 , L2.v)) )
        @test all( gety(filter(x -> x.x > 7.0 , e.env.v)) .== getx(filter(x -> x.x > 7.0 , L1.v)) )

        upper_env!(e, do_intersect = true)
        @test issorted(getx(e.env))
        @test gets(e)[1].x ≈ 5.0/0.7
        @test gets(e)[1].y ≈ 5.0/0.7
    end
    @testset "upper_env test 4" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        x2 = collect(range(-1,stop = 9,length = n))
        L1 = MLine(x1,x1)
        L2 = MLine(x2,5.0 .+ 0.3*x2)
        e = Envelope([L1,L2])
        upper_env!(e, do_intersect=true)
        @test issorted(getx(e.env))
        # @test isapprox(getx(e.env),sort(unique(vcat(x1,x2))), atol = 1e-15)
        @test gets(e)[1].x ≈ 5.0/0.7
        @test gets(e)[1].y ≈ 5.0/0.7
    end
    # poor test since we only do linear interpolation...
    @testset "log(x) vs 2(x-1)" begin
        x = collect(0.0:0.1:1.5)
        L1 = MLine(x,log.(x))
        L2 = MLine(x, 2 .* (x .- 1.0) )
        e = Envelope([L1,L2])
        upper_env!(e, do_intersect=true)
        @test issorted(getx(e.env))
        # @test gets(e)[1].x ≈ DCEGM.fzero(x-> 2*(x-1) - log(x),0.01,0.5)
        # @test gets(e)[1].y ≈ log( gets(e)[1].x )
        @test gets(e)[2].x == 1.0
        @test gets(e)[2].y == 0.0
    end
    @testset "y=x, y=2x, y=0.1x + 1.98" begin
        x1 = collect(-0.9:0.3:2.7)
        L1 = MLine(x1, x1)
        x2 = collect(0.0:0.1:1)
        L2 = MLine(x2, 2 .* x2, extrap = false)
        x3 = collect(1.0:0.45:2.9)
        L3 = MLine(x3, (0.1 .* x3) .+ 1.9, extrap = false)
        e = Envelope([L1,L2,L3])
        upper_env!(e)
        @test issorted(getx(e.env))
        @test all( gety(filter(x -> x.x < 0.0 , e.env.v)) .== gety(filter(x -> x.x < 0.0 , L1.v)) )
        @test all( gety(filter(x -> (x.x > 0.0) & (x.x < 1.0) , e.env.v)) .== gety(filter(x -> (x.x > 0.0) & (x.x < 1.0)  , L2.v)) )

        upper_env!(e,do_intersect = true)
        @test issorted(getx(e.env))
        @test all( gety(filter(x -> x.x < 0.0 , e.env.v)) .== gety(filter(x -> x.x < 0.0 , L1.v)) )
        @test all( gety(filter(x -> (x.x > 0.0) & (x.x < 1.0)  , e.env.v)) .== gety(filter(x -> (x.x > 0.0) & (x.x < 1.0)  , L2.v)) )
        # @test all( gety(filter(x -> (x.x >= 1.0)  , e.env.v)) .== gety(filter(x -> (x.x >= 1.0)  , L3.v)) )
        @test gets(e)[1] == Point(0.0,0.0)
        @test gets(e)[2] == Point(1.0,2.0)
        @test_broken gets(e)[3] == Point(1.9/0.9,1.9/0.9)   # something weird going on with intersect()



    end
    @testset "upper_env test: decreasing " begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        x2 = collect(range(-1,stop = 9,length = n))
        x = vcat(x1,x2)
        y = vcat(x1[end:-1:1],ones(n)*5)
        L = MLine(x,y)
        e = splitLine(L)
        upper_env!(e,do_intersect = true)
        @test issorted(getx(e.env))
        @test isapprox(getx(e.env), unique(sort(vcat(x1,x2))), atol = 1e-10)
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
    end
end
# @testset "Testing Envelopes over more MLines" begin
#     @testset "upper_env test 1" begin
#         n = 15
#         x1 = collect(range(0,stop = 10,length = n))
#         x2 = collect(range(-1,stop = 9,length = n))
#         x3 = collect(range(-0.1,stop = 10,length = n))
#         L1 = MLine(x1,x1)
#         L2 = MLine(x2,ones(n)*5)
#         L3 = MLine(x3,(x3.^2)/8)
#         e = Envelope([L1,L2,L3])
#         upper_env!(e)
#         @test issorted(getx(e.env))
#         @test getx(e.env) == sort(unique(vcat(x1,x2,x3,gets(e)[2].x)))
#         @test length(gets(e)) == 2
#         @test gets(e)[1].x ≈ 5.0
#         @test gets(e)[1].y ≈ 5.0
#         @test isapprox(gets(e)[2].x,8.0,rtol=0.01)
#         @test isapprox(gets(e)[2].y,8.0,rtol=0.01)
#     end
#     @testset "upper_env test 2" begin
#         n = 10
#         x1 = collect(range(1,stop = 10,length = n))
#         L1 = MLine(x1,x1.*(x1.<6))
#         L2 = MLine(x1,ones(n)*5 .* (x1.>4))
#         L3 = MLine(x1,x1 .- 3)
#         e = Envelope([L1,L2,L3])
#         upper_env!(e)
#         @test issorted(getx(e.env))
#         @test getx(e.env) == sort(unique(vcat(x1,[gets(e)[i].x for i in 1:2])))
#         @test length(gets(e)) == 2
#         @test gets(e)[1].x == 5.0
#         @test gets(e)[1].y == 5.0
#         @test gets(e)[2].x == 8.0
#         @test gets(e)[2].y == 5.0
#     end
# end

@testset "splitLine" begin
    @testset "test 1: simple" begin
        x = [1,2,3,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.2,1.8,2.1]
        L1 = MLine(x,y)
        e = splitLine(L1)
        @test isa(e,Envelope)
        @test length(getx(e.env))==1
        @test length(gety(e.env))==1
        @test length(gets(e))==0
        @test length(getr(e))==0
        @test length(e.L) == 3-1
        @test eltype(e.L)== MLine{Float64}

        @test extrema(getx(e.L[1])) == (1.0,3.0)
        @test extrema(getx(e.L[2])) == (1.5,2.9)
    end
    @testset "test 2: less simple" begin
        x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
        L1 = MLine(x,y)
        e = splitLine(L1)
        @test isa(e,Envelope)
        @test length(getx(e.env))==1
        @test length(gety(e.env))==1
        @test length(gets(e))==0
        @test length(e.L) == 3
        @test eltype(e.L)== MLine{Float64}

        @test extrema(getx(e.L[1].v)) == (1.0,3.0)
        @test extrema(getx(e.L[2].v)) == (1.5,3.0)
        @test extrema(getx(e.L[3].v)) == (1.5,2.9)
    end
    @testset "4 MLines - fails!" begin
        f1(x) = ones(length(x))
        f2(x) = 0.5x
        f3(x) = x .- 2
        f4(x) = 2x .- 8
        x1 = collect(range(0.9,stop = 1.9,length = 14))
        x2 = collect(range(1,stop = 7,length = 19))
        x3 = collect(range(1,stop = 7,length = 15))
        x4 = collect(range(1,stop = 8,length = 25))
        X = [x1...,x2...,x3...,x4...]
        L = MLine([x1...,x2...,x3...,x4...],vcat(f1(x1),f2(x2),f3(x3),f4(x4)))
        en = splitLine(L)

        @test isa(en,Envelope)
        @test length(getx(en.env))==1
        @test length(gety(en.env))==1
        @test length(gets(en))==0
        @test length(en.L) == 7-3

        upper_env!(en, do_intersect = true)
        @test issorted(getx(en.env))
        @test length(gets(en)) == 3
        xx = sort(vcat(X,vcat([gets(en)[i].x for i in 1:3]...)))
        @test isapprox(getx(en.env), unique(xx) , atol = 1e-10)
        yy = reshape(vcat([[f1(ix) f2(ix) f3(ix) f4(ix)] for ix in xx]...),length(xx),4)
        # @test isapprox(gety(en)[:],findmax(yy,2)[1][:])


    end

    @testset "4 MLines - passes!" begin
        f1(x) = ones(length(x))
        f2(x) = 0.5x
        f3(x) = x .-2
        f4(x) = 2x .-8
        x1 = collect(range(0.9,stop = 2.1,length = 14))
        x2 = collect(range(1,stop = 7,length = 19))
        x3 = collect(range(1,stop = 7,length = 15))
        x4 = collect(range(1,stop = 8,length = 25))
        X = [x1...,x2...,x3...,x4...]
        L = MLine([x1...,x2...,x3...,x4...],vcat(f1(x1),f2(x2),f3(x3),f4(x4)))
        en = splitLine(L)

        @test isa(en,Envelope)
        @test length(getx(en.env))==1
        @test length(gety(en.env))==1
        @test length(gets(en))==0
        @test length(en.L) == 7-3

        upper_env!(en, do_intersect = true)
        @test issorted(getx(en.env))
        @test isapprox(gets(en)[1] , Point(2.0,1.0))
        @test isapprox(gets(en)[2] , Point(4.0,2.0))
        @test isapprox(gets(en)[3] , Point(6.0,4.0))
        xx = unique(sort(vcat(X,vcat([gets(en)[i].x for i in 1:3]...))))
        @test isapprox(getx(en.env) , xx, atol = 1e-10)
        yy = reshape(vcat([[f1(ix) f2(ix) f3(ix) f4(ix)] for ix in xx]...),length(xx),4)
        @test isapprox(gety(en.env)[:],findmax(yy,dims=2)[1][:])


    end


end

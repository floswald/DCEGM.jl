

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
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2)))
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
    end
    @testset "upper_env test 2" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        L1 = MLine(x1,x1)
        L2 = MLine(x1,ones(n)*5)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test getx(e) == x1
        @test issorted(getx(e))
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
        @test getx(e) == x1
        @test issorted(getx(e))
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
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2)))
        @test gets(e)[1].x ≈ 5.0/0.7
        @test gets(e)[1].y ≈ 5.0/0.7
    end
    @testset "upper_env test: decreasing " begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        x2 = collect(range(-1,stop = 9,length = n))
        x = vcat(x1,x2)
        y = vcat(x1[end:-1:1],ones(n)*5)
        L = MLine(x,y)
        e = splitMLine(L)
        upper_env!(e)
        @test issorted(getx(e))
        @test_broken getx(e) == unique(sort(vcat(x1,x2)))
        @test gets(e)[1].x ≈ 0.0
        @test gets(e)[1].y ≈ 10.0
    end
end
@testset "Testing Envelopes over more MLines" begin 
    @testset "upper_env test 1" begin
        n = 15
        x1 = collect(range(0,stop = 10,length = n))
        x2 = collect(range(-1,stop = 9,length = n))
        x3 = collect(range(-0.1,stop = 10,length = n))
        L1 = MLine(x1,x1)
        L2 = MLine(x2,ones(n)*5)
        L3 = MLine(x3,(x3.^2)/8)
        e = Envelope([L1,L2,L3])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2,x3,gets(e)[2].x)))
        @test length(gets(e)) == 2
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
        @test isapprox(gets(e)[2].x,8.0,rtol=0.01)
        @test isapprox(gets(e)[2].y,8.0,rtol=0.01)
    end
    @testset "upper_env test 2" begin
        n = 10
        x1 = collect(range(1,stop = 10,length = n))
        L1 = MLine(x1,x1.*(x1.<6))
        L2 = MLine(x1,ones(n)*5 .* (x1.>4))
        L3 = MLine(x1,x1 .- 3)
        e = Envelope([L1,L2,L3])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,[gets(e)[i].x for i in 1:2])))
        @test length(gets(e)) == 2
        @test gets(e)[1].x == 5.0
        @test gets(e)[1].y == 5.0
        @test gets(e)[2].x == 8.0
        @test gets(e)[2].y == 5.0
    end
end

@testset "splitMLine" begin
    @testset "test 1: simple" begin
        x = [1,2,3,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.2,1.8,2.1]
        L1 = MLine(x,y)
        e = splitMLine(L1)
        @test isa(e,Envelope)
        @test length(getx(e))==1
        @test length(gety(e))==1
        @test length(gets(e))==0
        @test length(getr(e))==0
        @test length(e.L) == 3-1
        @test eltype(e.L)== MLine{Float64}

        @test extrema(e.L[1].x) == (1.0,3.0)
        @test extrema(e.L[2].x) == (1.5,2.9)
    end
    @testset "test 2: less simple" begin
        x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
        L1 = MLine(x,y)
        e = splitMLine(L1)
        @test isa(e,Envelope)
        @test length(getx(e))==1
        @test length(gety(e))==1
        @test length(gets(e))==0
        @test length(e.L) == 3
        @test eltype(e.L)== MLine{Float64}

        @test extrema(e.L[1].x) == (1.0,3.0)
        @test extrema(e.L[2].x) == (1.5,3.0)
        @test extrema(e.L[3].x) == (1.5,2.9)
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
        en = splitMLine(L)
 
        @test isa(en,Envelope)
        @test length(getx(en))==1
        @test length(gety(en))==1
        @test length(gets(en))==0
        @test length(en.L) == 7-3

        upper_env!(en)
        @test issorted(getx(en))
        @test length(gets(en)) == 3
        xx = sort(vcat(X,vcat([gets(en)[i].x for i in 1:3]...)))
        @test getx(en) == unique(xx)
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
        en = splitMLine(L)
 
        @test isa(en,Envelope)
        @test length(getx(en))==1
        @test length(gety(en))==1
        @test length(gets(en))==0
        @test length(en.L) == 7-3

        upper_env!(en)
        @test issorted(getx(en))
        @test gets(en) == [Point(2.0,1.0,i=1),Point(4.0,2.0,i=2),Point(6.0,4.0,i=3)]
        xx = unique(sort(vcat(X,vcat([gets(en)[i].x for i in 1:3]...))))
        @test getx(en) == xx
        yy = reshape(vcat([[f1(ix) f2(ix) f3(ix) f4(ix)] for ix in xx]...),length(xx),4)
        @test isapprox(gety(en)[:],findmax(yy,dims=2)[1][:])


    end


end




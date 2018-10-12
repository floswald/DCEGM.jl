

@testset "setup" begin
    x = collect(0:0.1:1)
    y = rand(11)
    m = Line(x,y)
    @test m.n == 11
    @test size(m) == (11,)
end

@testset "Interpolations" begin
    @testset "interpolate at log(1)" begin
        x = collect(0:0.1:1)
        y = log.(1 .+ x)
        m = Line(x,y)
        i = interp(m,[0.0])
        @test i[1] == 0.0
    end
end

@testset "getter and setters" begin
    x = [1,3,4,10]
    y = [6,0,8,1]
    l = Line(x,y)
    @test size(l) == (4,)

    @test l[1] == (1,6)
    @test isa(l[1:2],Line{Int64})
    @test l[1:2].x == [1,3]
    @test l[1:2].y == [6,0]

    setindex!(l,-1,-1,1)
    @test size(l) == (4,)
    @test l[1] == (-1,-1)

    setindex!(l,[10,11],[90,91],3:4)
    @test l[3:4].x == [10,11]
    @test l[3:4].y == [90,91]
    
end

@testset "Modifying Line methods" begin

    @testset "prepend" begin
        x = collect(0:0.1:1)
        y = log.(1 .+ x)
        L = Line(x,y)
        prepend!(L,18.0,-1.1)
        @test size(L)==(12,)
        @test L[1] == (18.0,-1.1)
        @test L[2] == (0.0,0.0)
    end
    
    @testset "append" begin
        x = collect(0:0.1:1)
        y = log.(1 .+ x)
        L = Line(x,y)
        append!(L,18.0,-1.1)
        @test size(L)==(12,)
        @test L[end] == (18.0,-1.1)
        @test L[1] == (0.0,0.0)
    end

    @testset "delete!" begin
        x = collect(0:0.1:1)
        y = rand(11)
        L = Line(x,y)
        delete!(L,9)
        @test size(L)==(10,)
        @test L[9] == (x[10],y[10])
    end
    @testset "insert!" begin
        x = collect(0:0.1:1)
        y = log.(1 .+ x)
        L = Line(x,y)
        insert!(L,1.1,3.3,9)
        @test size(L)==(12,)
        @test L[9] == (1.1,3.3)
        @test L[10] == (x[9],y[9])
        @test L[8] == (x[8],y[8])
    end
    @testset "splitat" begin
        x = collect(0:0.1:1)
        y = log.(1 .+ x)
        L = Line(x,y)
        n,o=splitat(L,1)
        @test size(n)==(1,)
        @test size(o)==(11,)

        n,o=splitat(L,2)
        @test size(n)==(2,)
        @test size(o)==(10,)
        @test n[2] == L[2]
        @test o[2] == L[3]
        @test o[1] == L[2]
        @test o[end] == L[end]

        n,o=splitat(L,2,false)
        @test size(n)==(2,)
        @test size(o)==(9,)
        @test n[2] == L[2]
        @test o[2] == L[4]
        @test o[1] == L[3]
    end
    @testset "sort!" begin
        x = collect(0:0.1:1)
        y = rand(11)
        L = Line(x,y)
        insert!(L,2.0,2.0,2)
        @test !(issorted(L.x))
        sort!(L)
        @test issorted(L.x)
    end
end



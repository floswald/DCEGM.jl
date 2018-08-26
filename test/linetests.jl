

@testset "Point Arithmetic" begin
    
    p = Point(1,1)
    @test eltype(p) == Int
    @test_throws MethodError Point(1,1.0) 

    p1 = Point(2,2)
    @test p + p1 == Point(3,3)

    p2 = Point(2.0,2.0)
    @test p + p1 == Point(3.0,3.0)

    @test 2*p2 == Point(4.0,4.0)
    @test p2 == Point(2.0,2.0)
    @test p2 != Point(2.0,1.9999999)

    arr = vcat(p2,Point(1.0,1.9999999))
    @test in(p2,arr)

end

@testset "Line construction from vectors" begin
    x = collect(0:0.1:1)
    y = rand(11)
    m = Line(x,y)
    @test m.n == 11
    @test size(m) == (11,)
    @test m.xvec == x
    @test m.yrange == extrema(y)
end

@testset "Line construction from vector of points" begin
    p1 = Point(1.0,2.3)
    p2 = Point(2.0,2.0)
    p3 = Point(3.0,5.2)
    m2 = Line([p1;p2;p3])
    @test length(m2) == 3
    @test size(m2) == (3,)
    @test m2.xrange == (1.0,3.0)
    @test m2.xvec == [1.0;2.0;3.0]
    @test m2.yrange == (2.0,5.2)
    @test eltype(m2) == Point{Float64}
    @test m2[1] == p1
    @test m2[3] == p3
    @test DCEGM.min_x(m2) == 1.0
    @test DCEGM.max_x(m2) == 3.0
    @test DCEGM.min_y(m2) == 2.0
    @test DCEGM.max_y(m2) == 5.2
end

# @testset "Interpolations" begin
#     @testset "interpolate at log(1)" begin
#         x = collect(0:0.1:1)
#         y = log.(1.+x)
#         m = Line(x,y)
#         i = interp(m,[0.0])
#         @test i[1] == 0.0
#     end
# end



@testset "Modifying Line methods" begin

    @testset "prepend" begin
        x = collect(0:0.1:1)
        y = log.(1.+x)
        L = Line(x,y)
        prepend!(L,[Point(18.0,-1.1)])
        @test size(L)==(12,)
        @test L[1] == Point(18.0,-1.1)
        @test L[2] == Point(0.0,0.0)
    end
    
    @testset "append" begin
        x = collect(0:0.1:1)
        y = log.(1.+x)
        L = Line(x,y)
        append!(L,[Point(18.0,-1.1)])
        @test size(L)==(12,)
        @test L[end] == Point(18.0,-1.1)
        @test L[1] == Point(0.0,0.0)
    end

    @testset "delete!" begin
        x = collect(0:0.1:1)
        y = rand(11)
        L = Line(x,y)
        delete!(L,9)
        @test size(L)==(10,)
        @test L[9] == Point(x[10],y[10])
        @test L.xvec == vcat(x[1:8],x[10:11])
    end
    @testset "insert!" begin
        x = collect(0:0.1:1)
        y = log.(1.+x)
        L = Line(x,y)
        insert!(L,Point(1.1,3.3),9)
        @test size(L)==(12,)
        @test L[9] == Point(1.1,3.3)
        @test L[10] == Point(x[9],y[9])
        @test L[8] == Point(x[8],y[8])
    end
    @testset "splitat" begin
        x = collect(0:0.1:1)
        y = log.(1.+x)
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
    @testset "sortx!" begin
        x = collect(0:0.1:1)
        y = rand(11)
        L = Line(x,y)
        insert!(L,Point(2.0,2.0),2)
        @test !(issorted(L))
        DCEGM.sortx!(L)
        @test issorted(L)
        @test issorted(L.xvec)
    end
end



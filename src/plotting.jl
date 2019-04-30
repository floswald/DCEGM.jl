
@recipe function f(m::Model2;id=1,iy=2)
    grid --> true
    xticks := true
    legend --> true
    for i in 1:(size(m.v)[3]-1)
        @series begin
            linetype --> :path 
            linewidth --> 1
            legend --> :bottomright
            seriescolor --> ColorGradient(:magma)
            m.v[id,iy,i].env.x,m.v[id,iy,i].env.y
        end
    end
end


struct tester
    m :: Matrix
end

@recipe function f(tt::tester)
    n,m = size(tt.m)
    for i in 1:n
        @series begin
            linetype --> :path
            series_annotations := ["$i" for i in 1:m]
            (1:m,tt.m[i,:])
        end
    end
end


@recipe function f(l::MLine;numerate=false,marker=false)
    # defaults
    grid --> true
    xticks := true
    legend := false
    @series begin
        linetype --> :path 
        linecolor --> :black
        linewidth --> 1
        if marker
            markershape --> :circle
            markerstrokecolor --> :black
            markercolor --> :white
            markersize --> 2
        end
        if numerate
            series_annotations := ["$i" for i in 1:l.n]
        end
        (getx(l),gety(l))
    end
end



# @recipe function f(x::Envelope; removed=false,num=false,marker=false)
@recipe function f(x::Envelope;numerate=false,removed=false,mrk=true)

    # defaults
    grid --> true
    xticks := true

    if !x.env_set
        legend --> :bottomright
    else
        legend --> false
    end

    println("marker = $mrk, title = $(get!(plotattributes,:title,""))")
    # println("marker = $mrk, numerate = $numerate, title = $(get!(plotattributes,:title,""))")

    # if line array exists, plot
    if length(x.L) > 0
        for l in x.L
            @series begin
                # subplot := 1
                seriestype := :path
                linewidth := 1
                if mrk
                    println("illegal")
                    markershape --> :circle
                    markerstrokecolor --> :black
                    markercolor --> :white
                    markersize := 3
                end
                if numerate
                    series_annotations := ["$i" for i in sortperm(getx(l))]
                end
                (getx(l),gety(l))
            end
        end
    end
    # plot envelope, if exists
    if x.env_set
        if removed
            for l in 1:length(x.L)
                ir = x.removed[l]
                if length(ir) > 0
                    @series begin
                        seriestype = :scatter
                        markershape := :rect
                        markersize := 3
                        markerstrokecolor := :black
                        markercolor := :white
                        markeralpha := 0.5
                        (getx(x.L[l])[ir],gety(x.L[l])[ir])
                    end
                end
            end
        end
        @series begin
            # subplot := 1
            seriestype := :path
            linecolor --> :red
            linewidth --> 2
            if mrk
                markershape := :circle
                markercolor := :white
                # markeralpha := 0.5
                markerstrokecolor := :black
                markersize := 3
            end
            if numerate
                series_annotations := ["$i" for i in 1:length(getx(x.env))]
            end
            (getx(x.env),gety(x.env))
        end
    end
end


function f3a()

    fs = Function[x->ones(length(x));x->0.5x;x->x .- 2;x->2x .- 8]
    xs = [
         collect(range(-1 , stop = 0.9 ,length= 6))   ,
         collect(range(1  , stop = 7   ,length= 19))  ,
         collect(range(2  , stop = 7   ,length= 15))  ,
         collect(range(4  , stop = 8   ,length= 25))]
    ls = [MLine(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    e = Envelope(ls)
    return e
end
function f3b()

    fs = Function[x->ones(length(x));x->0.5x;x->x .- 2;x->2x .- 8]
    xs = [
         collect(range(-1 , stop = 3.1 , length = 6))   ,
         collect(range(1  , stop = 7   , length = 19))  ,
         collect(range(2  , stop = 7   , length = 15))  ,
         collect(range(4  , stop = 8   , length = 25))]
    ls = [MLine(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    e = Envelope(ls)
    return e
end
function f3c()

    fs = Function[x->ones(length(x));x -> 0.5x;x->x .- 2;x -> 2x .- 8.5]
    xs = [
         collect(range(-1 , stop = 3.1 , length= 6))   ,
         collect(range(1  , stop = 7   , length= 19))  ,
         collect(range(2  , stop = 7   , length= 15))  ,
         collect(range(4  , stop = 8   , length= 25))]
    ls = [MLine(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    e = Envelope(ls)
    return e
end

function tplot1()
    n = 15
    x1 = collect(range(0 , stop = 10,length = n))
    x2 = collect(range(-1, stop = 9 ,length = n))
    L1 = MLine(x1,x1)
    L2 = MLine(x2,ones(n)*5)
    e = Envelope([L1,L2])
    plot(e)
end

function tplot_intersect()
    n = 15
    x1 = collect(range(0 , stop = 10,length = n))
    x2 = collect(range(-0.5, stop = 9 ,length = n))
    L1 = MLine(x1,x1)
    L2 = MLine(x2,ones(n)*4.6)
    e = Envelope([L1,L2])
    p1 = plot(e)
    upper_env!(e)
    p2 = plot(e,title = "do_intersect = false",mrk = false)
    upper_env!(e,do_intersect = true)
    p3 = plot(e,title = "do_intersect = true",mrk=false)
    plot(p1,p2,p3,layout = (1,3))
end

function tplot2()
    n = 15
    x1 = collect(range(0  , stop=10 ,length= n))
    x2 = collect(range(-1 , stop=9  ,length= n))
    L1 = MLine(x1,x1)
    L2 = MLine(x2,ones(n)*5)
    en = Envelope([L1,L2])
    p1 = plot(en)

    upper_env!(en)
    removed!(en)
    p2 = plot(en,removed=true)
    plot(p1,p2)
end
function tplot3a()

    en = f3a()

    p1 = plot(en,title = "non-overlapping grids")

    upper_env!(en)
    removed!(en)
    p2 = plot(en)
    p2 = plot(en,removed=true,title = "Envelope and removed points")
    plot(p1,p2)

end

function tplot3b()

    en = f3b()

    p1 = plot(en,title = "overlapping grids")

    upper_env!(en)
    removed!(en)
    p2 = plot(en,removed=true)

    plot(p1,p2)

end
function tplot3c()

    en = f3c()

    p1 = plot(en)

    upper_env!(en)
    removed!(en)
    p2 = plot(en,removed=true)

    plot(p1,p2)
    gui()
    return en

end

function tplot4()

    fs = Function[x->ones(length(x));x->0.5x;x->x .- 2;x->2x .- 8]
    xs = [
         collect(range(0.1 , stop =  1.5 ,length = 5))   ,
         collect(range(1   , stop =  7   ,length = 19))  ,
         collect(range(2   , stop =  7   ,length = 15))  ,
         collect(range(4   , stop =  8   ,length = 25))]
    ls = [MLine(i[1],i[2](i[1])) for i in zip(xs,fs)]
    e = Envelope(ls)

    p1 = plot(e)

    upper_env!(e)
    removed!(e)

    p2 = plot(e,removed=true)

    plot(p1,p2)
end

function tplot5()

    x1 = collect(-0.9:0.3:2.7)
    L1 = MLine(x1, x1)
    x2 = collect(0.0:0.1:1)
    L2 = MLine(x2, 2 .* x2, extrap = false)
    x3 = collect(1.0:0.45:2.9)
    L3 = MLine(x3, (0.1 .* x3) .+ 1.9, extrap = false)
    e = Envelope([L1,L2,L3])
    p1 = plot(e,title = "set extrap=false")
    upper_env!(e)
    removed!(e)
    p2 = plot(e,removed=true, title = "correct envelope")
    println("env = $(e.env.v)")
    plot(p1,p2)

end

function splitf()
    x = [1,2,3,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.2,1.8,2.1]
    L = MLine(x,y)
    p1 = plot(L,title="original",numerate=true)
    e = splitLine(L)
    p2 = plot(e,title="split MLine",numerate=true)
    plot(p1,p2)
end

function splitf2()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = MLine(x,y)
    p1 = plot(L,title="original",numerate=true)
    e = splitLine(L)
    p2 = plot(e,title="split MLine",numerate=true)
    plot(p1,p2)
end

function splitf3()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = MLine(x,y, extrap = false)
    e = splitLine(L)
    upper_env!(e)
    p1 = plot(e,title="MLine(x,y,extrap=false)")
    removed!(e)
    p2 = plot(e,title="with removed points",removed=true)
    plot(p1,p2)
end



function allplots()
    p = joinpath(dirname(@__FILE__),"..","images")
    tplot1()
    savefig(joinpath(p,"tplot1.png"))
    tplot2()
    savefig(joinpath(p,"tplot2.png"))
    tplot3a()
    savefig(joinpath(p,"tplot3a.png"))
    tplot3b()
    savefig(joinpath(p,"tplot3b.png"))
    tplot3c()
    savefig(joinpath(p,"tplot3c.png"))
    tplot5()
    savefig(joinpath(p,"tplot5.png"))
    tplot_intersect()
    savefig(joinpath(p,"tplot_intersect.png"))
    splitf()
    savefig(joinpath(p,"split.png"))
    splitf2()
    savefig(joinpath(p,"split2.png"))
    splitf3()
    savefig(joinpath(p,"split3.png"))
end




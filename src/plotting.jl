
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


@recipe function f(l::Line;numerate=false,marker=false)
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
            series_annotations := ["$i" for i in 1:length(l.x)]
        end
        (l.x,l.y)
    end
end



# @recipe function f(x::Envelope; removed=false,num=false,marker=false)
@recipe function f(x::Envelope;numerate=false,removed=false,marker=true)

    # defaults
    grid --> true
    xticks := true

    if !x.env_set
        legend --> :bottomright
    else
        legend --> false
    end

    # if line array exists, plot
    if length(x.L) > 0
        for l in x.L
            @series begin
                # subplot := 1
                linetype := :path
                linewidth := 1
                if marker
                    markershape --> :circle
                    markerstrokecolor --> :black
                    markercolor --> :white
                    markersize := 3
                end
                if numerate
                    series_annotations := ["$i" for i in sortperm(l.x)]
                end
                (getx(l),gety(l))
            end
        end
    end
    # plot envelope, if exists
    if x.env_set
        @series begin
            # subplot := 1
            linetype := :line 
            linecolor --> :red
            linewidth --> 2
            if marker
                markershape := :circle
                markercolor := :white
                # markeralpha := 0.5
                markerstrokecolor := :black
                markersize := 3
            end
            if numerate
                series_annotations := ["$i" for i in 1:length(getx(x))]
            end
            (getx(x.env),gety(x.env))
        end
        if removed
            for l in 1:length(x.L)
                @series begin
                    seriestype = :scatter
                    markershape := :rect
                    markersize := 3
                    markerstrokecolor := :black
                    markercolor := :white
                    markeralpha := 0.5
                    (x.L[l].x[x.removed[l]], x.L[l].y[x.removed[l]])
                end
            end
        end
    end
end

function tplot1()
    n = 15
    x1 = collect(linspace(0,10,n))
    x2 = collect(linspace(-1,9,n))
    L1 = Line(x1,x1)
    L2 = Line(x2,ones(n)*5)
    e = Envelope([L1,L2])
    plot(e)
end

function tplot2()
    n = 15
    x1 = collect(linspace(0,10,n))
    x2 = collect(linspace(-1,9,n))
    L1 = Line(x1,x1)
    L2 = Line(x2,ones(n)*5)
    en = Envelope([L1,L2])
    p1 = plot(en)

    upper_env!(en)
    p2 = plot(en,removed=true)
    plot(p1,p2)
end

function f3a()

    fs = Function[x->ones(length(x));x->0.5x;x->x-2;x->2x-8]
    xs = [
         collect(linspace(-1,0.9,6)),
         collect(linspace(1,7,19)),
         collect(linspace(2,7,15)),
         collect(linspace(4,8,25))]
    ls = [Line(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    e = Envelope(ls)
    return e
end
function f3b()

    fs = Function[x->ones(length(x));x->0.5x;x->x-2;x->2x-8]
    xs = [
         collect(linspace(-1,3.1,6)),
         collect(linspace(1,7,19)),
         collect(linspace(2,7,15)),
         collect(linspace(4,8,25))]
    ls = [Line(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    e = Envelope(ls)
    return e
end
function f3c()

    fs = Function[x->ones(length(x));x->0.5x;x->x-2;x->2x-8.5]
    xs = [
         collect(linspace(-1,3.1,6)),
         collect(linspace(1,7,19)),
         collect(linspace(2,7,15)),
         collect(linspace(4,8,25))]
    ls = [Line(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    e = Envelope(ls)
    return e
end

function tplot3a()

    en = f3a()

    p1 = plot(en)

    upper_env!(en)
    p2 = plot(en,removed=true)
    plot(p1,p2)

end

function tplot3b()

    en = f3b()

    p1 = plot(en)

    upper_env!(en)
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

    fs = Function[x->ones(length(x));x->0.5x;x->x-2;x->2x-8]
    xs = [
         collect(linspace(0.1,1.5,5)),
         collect(linspace(1,7,19)),
         collect(linspace(2,7,15)),
         collect(linspace(4,8,25))]
    ls = [Line(i[1],i[2](i[1])) for i in zip(xs,fs)]
    e = Envelope(ls)

    p1 = plot(e)

    upper_env!(e)

    p2 = plot(e,removed=true)

    plot(p1,p2)
end

function splitf()
    x = [1,2,3,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.2,1.8,2.1]
    L = Line(x,y)
    p1 = plot(L,title="original",numerate=true)
    e = splitLine(L)
    p2 = plot(e,title="split Line",numerate=true)
    plot(p1,p2)
end

function splitf2()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = Line(x,y)
    p1 = plot(L,title="original",numerate=true)
    e = splitLine(L)
    p2 = plot(e,title="split Line",numerate=true)
    plot(p1,p2)
end

function splitf3()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = Line(x,y)
    e = splitLine(L)
    upper_env!(e)
    p1 = plot(e,title="upper enveloped")
    removed!(e)
    p2 = plot(e,title="with removed points",removed=true)
    plot(p1,p2)
end

function allplots()
    p = joinpath(dirname(@__FILE__),"..","images")
    splitf2()
    savefig(joinpath(p,"split1.png"))
    splitf3()
    savefig(joinpath(p,"split2.png"))
end




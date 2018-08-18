
@recipe function f(l::Line;numerate=false,m=false)
    # defaults
    grid --> true
    xticks := true
    legend := false
    @series begin
        linetype --> :path 
        linecolor --> :black
        linewidth --> 1
        if m
            markershape --> :circle
            markerstrokecolor --> :black
            markercolor --> :white
            markersize --> 2
            if numerate
                series_annotations := ["$i" for i in 1:length(l.x)]
            end
        end
        (l.x,l.y)
    end
end



@recipe function f(x::Envelope; removed=false,numerate=false,m=false)

    # defaults
    grid --> true
    xticks := true
    legend := false

    # if line array exists, plot
    if length(x.L) > 0
        for l in x.L
            @series begin
                # subplot := 1
                linetype := :line 
                linecolor := :black
                linewidth := 1
                if m
                    markershape := :circle
                    markerstrokecolor := :black
                    markercolor := :white
                    markersize := 2
                    if numerate
                        series_annotations := ["$i" for i in 1:length(l.x)]
                    end
                end
                (l.x,l.y)
            end
        end
    end
    # plot envelope, if exists
    if x.env_set
        @series begin
            # subplot := 1
            linetype := :line 
            linecolor --> :red
            linewidth --> 4
            if m
                markershape := :circle
                markercolor := :white
                # markeralpha := 0.5
                markerstrokecolor := :black
                markersize := 3
                if numerate
                    series_annotations := ["$i" for i in 1:length(getx(x))]
                end
            end
            (getx(x),gety(x))
        end
        if removed
            for ir in x.removed
                if length(ir) > 0
                    @series begin
                        seriestype = :scatter
                        markershape := :rect
                        markersize := 3
                        markerstrokecolor := :black
                        markercolor := :white
                        markeralpha := 0.5
                        [ir[i].x for i in 1:length(ir)],[ir[i].y for i in 1:length(ir)]
                    end
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
    p2 = plot(en,removed=true)

    plot(p1,p2)

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




function v_analytic(m::Model,p::Param,id,it)
    vf = m.v[id,it]
    c = m.c[id,it]
    cons = scaleGrid(p.cfloor_plot,gety(c.env)[2],p.k,logorder = 1)
    cash = scaleGrid(p.a_low,getx(vf.env)[2],p.k,logorder = 1)
    # deleteat!(cons,length(cons))
    # deleteat!(cash,length(cash))
    pts = convert(Point,cash,vfun(id,it,cons,cash,vf,p))
    vcat(pts, vf.env.v[2:end])  # connect at second point
end

function v_analytic(m::Model,p::Param,id,iy,it)
    vf = m.v[id,iy,it]
    c = m.c[id,iy,it]
    cons = scaleGrid(p.cfloor_plot,gety(c.env)[2],p.k,logorder = 1)
    cash = scaleGrid(p.a_low,getx(vf.env)[2],p.k,logorder = 1)
    # deleteat!(cons,length(cons))
    # deleteat!(cash,length(cash))
    pts = convert(Point,cash,vfun(id,it,cons,cash,vf,p))
    vcat(pts, vf.env.v[2:end])  # connect at second point
    # vf.env.v[2:end]  # connect at second point
end

function plot_s(s::Simulation)

    # inc, cons, w
    py = plot(s.inc',leg = false, title = "income")
    pc = plot(s.cons',leg = false, title = "consumption")
    i_retires = findall(s.ret_age .> 0)
    scatter!(pc, s.ret_age[i_retires], [s.cons[i,s.ret_age[i]] for i in i_retires], m = (:rect, 2, 0.6, :white))
    pw0 = plot(s.w0',leg = false, title = "w0")
    scatter!(pw0, s.ret_age[i_retires], [s.w0[i,s.ret_age[i]] for i in i_retires], m = (:rect, 2, 0.6, :white))

    ppr = plot(s.prob_work',leg = false, title = "p(work)")
    scatter!(ppr, s.ret_age[i_retires], [s.prob_work[i,s.ret_age[i]] for i in i_retires], m = (:rect, 2, 0.6, :white))

    pw1 = plot(s.w1',leg = false, title = "w1")
    pempty = plot(legend=false,grid=false,foreground_color_subplot=:white)
    asize = 10
    annotate!(pempty, [(0.3,1,Plots.text("gamma = $(s.p.gamma)",     :left, asize)),
                       (0.3,0.9,Plots.text("R = $(s.p.R)",     :left, asize)),
                       (0.3,0.8,Plots.text("beta = $(round(s.p.beta,digits=2))", :left, asize)),
                       (0.3,0.7,Plots.text("alpha = $(s.p.alpha)",   :left, asize)),
                       (0.3,0.6,Plots.text("sigma = $(s.p.sigma)",   :left, asize)),
                       (0.3,0.5,Plots.text("lambda = $(s.p.lambda)", :left, asize)),
                       (0.3,0.4,Plots.text("rho = $(s.p.ρ)",         :left, asize))])
    plot(py,pc,pw0,ppr,pw1,
          pempty ,
          layout = (2,3))
end

@recipe function f(m::GModel,p::Param;id=1,iy=nothing,it=nothing)
    grid --> true
    xticks := true
    legend --> false
    # cg = cgrad(:inferno)
    # c1 = colorant"red"
    # c2 = colorant"blue"
    xrange = zeros(2)
    xa = extrema(m.avec)
    xrange[1] = xa[1] .- diff(vcat(xa...))[1] .* 0.01
    xrange[2] = xa[2]

    nT = size(m.v)[3]
    # cols = range(c1,stop=c2,length=nT)

    layout := grid(1,2)

    if isnothing(iy) & !isnothing(it)
        legend := true
        title --> ["value period $it" "consumption period $it"]
        for jy in 1:p.ny
            vt = v_analytic(m,p,id,jy,it)
            @series begin
                linetype --> :path
                linewidth --> 1
                legend := :bottomright
                label := "y$jy"
                # seriescolor --> cols[i]
                subplot := 1  # value function
                yguide := "value"
                xguide := "Cash on Hand M"
                xlim := xrange
                # ylim --> (-15,15)
                getx(vt),gety(vt)
                # getx(m.v[id,iy,it].env),gety(m.v[id,iy,it].env)
            end
        end
        for jy in 1:p.ny
            @series begin
                linetype --> :path
                linewidth --> 1
                legend --> :bottomright
                # seriescolor --> cols[i]
                subplot := 2  # 
                label := "y$jy"
                xlim := xa
                ylim --> xa
                yguide := "consumption"
                xguide := "Cash on Hand M"
                # aspect_ratio := aspect
                getx(m.c[id,jy,it].env),gety(m.c[id,jy,it].env)
            end
        end
    elseif isnothing(it) & !isnothing(iy)
        for i in 1:nT
            vt = v_analytic(m,p,id,iy,i)
            @series begin
                linetype --> :path
                linewidth --> 1
                legend --> :bottomright
                # seriescolor --> cols[i]
                # label := lab
                xlim := xrange
                # ylim --> (-15,15)
                subplot := 1  # value function
                yguide := "value"
                xguide := "Cash on Hand M"
                getx(vt),gety(vt)
            end
            @series begin
                linetype --> :path
                linewidth --> 1
                legend --> :bottomright
                # seriescolor --> cols[i]
                subplot := 2  # 
                # label := lab
                xlim := xa
                ylim --> xa
                yguide := "consumption"
                xguide := "Cash on Hand M"
                # aspect_ratio := aspect
                getx(m.c[id,iy,i].env),gety(m.c[id,iy,i].env)
            end
        end
    else
        println("you need to either give it or iy. not both. not none.")
    end
    # elseif isnothing(it) & isnothing(iy)
    #     title --> ["value period $it" "value period $it"]
    #     vt = v_analytic(m,p,id,iy,it)
    #     @series begin
    #         linetype --> :path
    #         linewidth --> 1
    #         legend --> :bottomright
    #         # seriescolor --> cols[i]
    #         subplot := 1  # value function
    #         yguide := "value"
    #         xguide := "Cash on Hand M"
    #         xlim := xrange
    #         # ylim --> (-15,15)
    #         getx(vt),gety(vt)
    #         # getx(m.v[id,iy,it].env),gety(m.v[id,iy,it].env)
    #     end
    #     @series begin
    #         linetype --> :path
    #         linewidth --> 1
    #         legend --> :bottomright
    #         # seriescolor --> cols[i]
    #         subplot := 2  # 
    #         xlim := xa
    #         ylim --> xa
    #         yguide := "consumption"
    #         xguide := "Cash on Hand M"
    #         # aspect_ratio := aspect
    #         getx(m.c[id,iy,it].env),gety(m.c[id,iy,it].env)
    #     end
    # end

end

@recipe function f(m::FModel,p::Param;id=1)
    grid --> true
    xticks := true
    legend --> true
    cg = cgrad(:inferno)
    # c1 = colorant"red"
    # c2 = colorant"blue"
    # alow,ahi = extrema(m.avec)
    # aspect = (ahi-alow)/(ahi - 0.0)

    nT = size(m.v)[2]
    # cols = range(c1,stop=c2,length=nT)
    cols = cg[range(0.0,stop=1.0,length=nT)]
    xrange = zeros(2)
    xa = extrema(m.avec)
    xrange[1] = xa[1] .- diff(vcat(xa...))[1] .* 0.01
    xrange[2] = xa[2]

    layout := grid(1,2)
    for i in 1:nT
        vt = v_analytic(m,p,id,i)
        lab = ((i==1)|(i==nT)) ? "$i" : ""
        @series begin
            linetype --> :path
            linewidth --> 1
            legend --> :bottomright
            seriescolor --> cols[i]
            label := lab
            subplot := 1  # value function
            yguide := "value"
            xguide := "Cash on Hand M"
            xlim := xrange
            ylim := (-15,15)
            getx(vt),gety(vt)
        end
        @series begin
            linetype --> :path
            linewidth --> 1
            legend --> :bottomright
            seriescolor --> cols[i]
            # subplot := 2  # 
            label := lab
            # xlim := (alow,ahi)
            # ylim := (0,ahi)
            ylim := extrema(m.avec)
            xlim := extrema(m.avec)
            yguide := "consumption"
            xguide := "Cash on Hand M"
            # aspect_ratio := aspect
            getx(m.c[id,i].env),gety(m.c[id,i].env)
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


@recipe function f(l::Vector{Point{T}};numerate=false,marker=false) where T
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
            markersize --> 1
        end
        if numerate
            series_annotations := ["$i" for i in 1:length(l.v)]
        end
        (getx(l),gety(l))
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
            markersize --> 1
        end
        if numerate
            series_annotations := ["$i" for i in 1:length(l.v)]
        end
        (getx(l),gety(l))
    end
end


@recipe function f(L::Vector{MLine{T}};numerate=false,mrk=true) where T
    for l in L
        @series begin
            # subplot := 1
            seriestype := :path
            linewidth := 1
            if mrk
                # println("illegal")
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
#
@recipe function f(x::Envelope;numerate=false,mrk=true)

    # defaults
    grid --> true
    xticks := true

    any_isec = length(x.isects) > 0

    # println("marker = $mrk, title = $(get!(plotattributes,:title,""))")
    # println("marker = $mrk, numerate = $numerate, title = $(get!(plotattributes,:title,""))")

    # actual envelope
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
    if any_isec && mrk
        @series begin
            seriestype := :scatter
            markershape := :star6
            markerstrokecolor := :black
            markercolor := :yellow
            markersize := 5
            (getx(x.isects),gety(x.isects))
        end
    end
end
#
@recipe function f(x::Envelope{T},L::Vector{MLine{T}};numerate=false,removed=false,mrk=true) where T

    # defaults
    grid --> true
    xticks := true
    legend --> false

    any_isec = length(x.isects) > 0

    # println("marker = $mrk, title = $(get!(plotattributes,:title,""))")
    # println("marker = $mrk, numerate = $numerate, title = $(get!(plotattributes,:title,""))")

    # lien arrya
    for l in L
        @series begin
            # subplot := 1
            seriestype := :path
            linewidth := 1
            if mrk
                # println("illegal")
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

    # actual envelope
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
    if any_isec && mrk
        @series begin
            seriestype := :scatter
            markershape := :star6
            markerstrokecolor := :black
            markercolor := :yellow
            markersize := 5
            (getx(x.isects),gety(x.isects))
        end
    end
end

# plot env and initial mline next to each other, showing new points and removed ones
@recipe function f(x::Envelope{T},L::MLine{T};numerate=false) where T

    # defaults
    grid --> true
    xticks := true
    legend --> false
    layout := (1,2)
    title --> ["Backwards-Bending" "Envelope"]

    any_isec = length(x.isects) > 0
    any_rmv = length(x.removed) > 0

    # println("marker = $mrk, title = $(get!(plotattributes,:title,""))")
    # println("marker = $mrk, numerate = $numerate, title = $(get!(plotattributes,:title,""))")

    # dirty line
    @series begin
        subplot := 1
        seriestype := :path
        linewidth := 1
        markershape --> :circle
        markerstrokecolor --> :black
        markercolor --> :white
        markersize := 3
        if numerate
            series_annotations := ["$i" for i in 1:length(getx(L))]
        end
        (getx(L),gety(L))
    end
    if any_rmv
        @series begin
            subplot := 1
            seriestype = :scatter
            markershape := :rect
            markersize := 3
            markerstrokecolor := :black
            markercolor := :white
            markeralpha := 0.5
            (getx(L)[x.removed],gety(L)[x.removed])
        end
    end

    # actual envelope
    @series begin
        subplot := 2
        seriestype := :path
        linecolor --> :red
        linewidth --> 1
        markershape --> :circle
        markerstrokecolor --> :black
        markercolor --> :white
        markersize := 3
        (getx(x.env),gety(x.env))
    end
    if any_isec
        @series begin
            subplot := 2
            seriestype := :scatter
            markershape := :star6
            markerstrokecolor := :black
            markercolor := :yellow
            markersize := 5
            (getx(x.isects),gety(x.isects))
        end
    end

end

#
# # @recipe function f(x::Envelope; removed=false,num=false,marker=false)
# @recipe function f(x::Envelope;numerate=false,removed=false,mrk=true)
#
#     # defaults
#     grid --> true
#     xticks := true
#
#     if !x.env_clean
#         legend --> :bottomright
#     else
#         legend --> false
#     end
#     any_isec = length(x.isects) > 0
#
#     # println("marker = $mrk, title = $(get!(plotattributes,:title,""))")
#     # println("marker = $mrk, numerate = $numerate, title = $(get!(plotattributes,:title,""))")
#
#     # if line array exists, plot
#     if length(x.L) > 0
#         for l in x.L
#             @series begin
#                 # subplot := 1
#                 seriestype := :path
#                 linewidth := 1
#                 if mrk
#                     # println("illegal")
#                     markershape --> :circle
#                     markerstrokecolor --> :black
#                     markercolor --> :white
#                     markersize := 3
#                 end
#                 if numerate
#                     series_annotations := ["$i" for i in sortperm(getx(l))]
#                 end
#                 (getx(l),gety(l))
#             end
#         end
#     end
#     # plot envelope, if exists
#     if x.env_clean
#         if removed
#             for l in 1:length(x.L)
#                 ir = x.removed[l]
#                 if length(ir) > 0
#                     @series begin
#                         seriestype = :scatter
#                         markershape := :rect
#                         markersize := 3
#                         markerstrokecolor := :black
#                         markercolor := :white
#                         markeralpha := 0.5
#                         (getx(x.L[l])[ir],gety(x.L[l])[ir])
#                     end
#                 end
#             end
#         end
#         @series begin
#             # subplot := 1
#             seriestype := :path
#             linecolor --> :red
#             linewidth --> 2
#             if mrk
#                 markershape := :circle
#                 markercolor := :white
#                 # markeralpha := 0.5
#                 markerstrokecolor := :black
#                 markersize := 3
#             end
#             if numerate
#                 series_annotations := ["$i" for i in 1:length(getx(x.env))]
#             end
#             (getx(x.env),gety(x.env))
#         end
#         if any_isec && mrk
#             @series begin
#                 seriestype := :scatter
#                 markershape := :star6
#                 markerstrokecolor := :black
#                 markercolor := :yellow
#                 markersize := 5
#                 (getx(x.isects),gety(x.isects))
#             end
#         end
#     end
# end




function f3a()

    fs = Function[x->ones(length(x));x->0.5x;x->x .- 2;x->2x .- 8]
    xs = [
         collect(range(-1 , stop = 0.9 ,length= 6))   ,
         collect(range(1  , stop = 7   ,length= 19))  ,
         collect(range(2  , stop = 7   ,length= 15))  ,
         collect(range(4  , stop = 8   ,length= 25))]
    ls = [MLine(i[1],i[2](i[1])) for i in zip(xs,fs)]
    return ls
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
    return ls
end
function f3c()

    fs = Function[x->ones(length(x));x -> 0.5x;x->x .- 2;x -> 2x .- 8.5]
    xs = [
         collect(range(-1 , stop = 3.1 , length= 6))   ,
         collect(range(1  , stop = 7   , length= 4))  ,
         collect(range(2  , stop = 7   , length= 15))  ,
         collect(range(4  , stop = 8   , length= 25))]
    ls = [MLine(i[1],i[2](i[1])) for i in zip(xs,fs)]
    # create an envelope
    return ls
end

function tplot1()
    n = 15
    x1 = collect(range(0 , stop = 10,length = n))
    x2 = collect(range(-1, stop = 9 ,length = n))
    L1 = MLine(x1,x1)
    L2 = MLine(x2,ones(n)*5)
    e = [L1,L2]
    x1,x2,e #,plot(e)
end

function tplot2()
    n = 15
    x1 = collect(range(0  , stop=10 ,length= n))
    x2 = collect(range(-1 , stop=9  ,length= n))
    L1 = MLine(x1,x1)
    L2 = MLine(x2,ones(n)*5)
    LL = [L1,L2]

    # a = splitLine(LL)
    e = upper_env(LL)
    p2 = plot(e,LL)
    savefig(p2,joinpath(@__DIR__,"..","images","tplot2.png"))
    p2
end
function tplot3a()

    en = f3a()

    p1 = plot(en,title = "non-overlapping grids")

    e = upper_env(en)
    p2 = plot(e,en,title = "Envelope")
    plot(p1,p2)

end

function tplot3b()

    en = f3b()

    p1 = plot(en,title = "overlapping grids")

    e = upper_env(en)
    p2 = plot(en)

    plot(p1,p2)

end
function tplot3c()

    en = f3c()

    p1 = plot(en)

    e = upper_env(en)
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
    e = DCEGM.upper_env(ls)

    p1 = plot(ls)
    p2 = plot(e,ls,title="envelope")
    p = plot(p1,p2)
    savefig(p,joinpath(@__DIR__,"..","images","tplot4.png"))
    p
end

function tplot5()

    x1 = collect(-0.9:0.3:2.7)
    L1 = MLine(x1, x1)
    x2 = collect(0.0:0.1:1)
    L2 = MLine(x2, 2 .* x2)
    x3 = collect(1.0:0.45:2.9)
    L3 = MLine(x3, (0.1 .* x3) .+ 1.9 )
    ls = [L1,L2,L3]
    e = DCEGM.upper_env(ls)
    p1 = plot(ls)
    p2 = plot(e,ls,title="envelope")

    p = plot(p1,p2)
    savefig(p,joinpath(@__DIR__,"..","images","tplot5.png"))
    p

end

function splitf()
    x = [1,2,3,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.2,1.8,2.1]
    L = MLine(x,y)
    p1 = plot(L,title="original",numerate=true)
    e = splitLine(L)
    p2 = plot(e,title="split MLine",numerate=true)
    p = plot(p1,p2)
    savefig(p,joinpath(@__DIR__,"..","images","split.png"))
    p
end

function splitf2()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = MLine(x,y)
    p1 = plot(L,title="original",numerate=true)
    e = splitLine(L)
    p2 = plot(e,title="split MLine",numerate=true)
    p = plot(p1,p2)
    savefig(p,joinpath(@__DIR__,"..","images","split2.png"))
    p
end

function splitf3()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = MLine(x,y)
    e = DCEGM.secondary_envelope(L)
    p2 = plot(e,L,numerate=true)
    savefig(p2,joinpath(@__DIR__,"..","images","split3.png"))
    p2
end

function splitf3a()
    x = [1,2,3,2.9,2.5,1.9,1.8,1.5,2.1,2.9]
    y = [1,1.5,1.7,1.6,1.55,1.4,1.3,1.2,1.8,2.1]
    L = MLine(x,y)
    e = secondary_envelope(L)
    plot(e,L)

end

function test_upper_env_dec()
    n = 15
    x1 = collect(range(0,stop = 10,length = n))
    x2 = collect(range(-1,stop = 9,length = n))
    x = vcat(x1,x2)
    y = vcat(x1[end:-1:1],ones(n)*5)
    L = MLine(x,y)
    e = splitLine(L)
    p1 = plot(L,numerate=true,title = "initial line")
    p2 = plot(e,marker=true,title = "split line")
    e = upper_env(e)
    p3 = plot(e,marker=true,title = "Envelope w new pts")
    p4 = plot(p1,p2,p3,layout = (1,3),size=(1200,500))
    savefig(p4,joinpath(@__DIR__,"..","images","descending_upper.png"))
    L,e,x1,x2,p4
end

function split_test()
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
    en,L
end

function demo(;n = 10,k = 10)
    a = Array{MLine{Float64}}(undef,10)
    for i in 1:2:n
        a[i] = MLine(sort([0.0; 10*rand(k); 10]), [0.0; rand(k+1).*collect(range(0.5,stop=10,length=k+1))])
        a[i+1] = MLine([0.0; 10],[10.0-i ; -(10-i)*1.5/i+10-i] )
    end
    p1 = plot(a,legend=false,title = "$n lines")
    e = upper_env(a)
    p2 = plot(e,a,title = "envelope")
    p =plot(p1,p2)
    savefig(p,joinpath(@__DIR__,"..","images","demo.png"))
    p
end

function demo2(;n = 10,k = 10)
    a = Array{MLine{Float64}}(undef,10)
    slope = 1.5
    for i in 1:2:n
        a[i] = MLine(sort([0.0; 10*rand(k); 10]), [0.0; rand(k+1).*collect(range(0.5,stop=10,length=k+1))])
        a[i+1] = MLine([0.0; 10],[10.0-i ; -(10-i)*slope/i+10-i] )
    end
    # add a particularly tricky line
    tricky = [a[2].v[1]]
    push!(tricky,DCEGM.Point(a[2].v[2].x,a[2].v[2].y + slope))
    push!(a, MLine(tricky))

    p1 = plot(a,legend=false,title = "$n lines")
    e = upper_env(a)
    p2 = plot(e,a,title = "first intersection on grid")
    p = plot(p1,p2)
    savefig(p,joinpath(@__DIR__,"..","images","demo2.png"))
    p,e
end



function allplots()
    tplot1()
    tplot2()
    # savefig(joinpath(p,"tplot3a.png"))
    # tplot3b()
    # savefig(joinpath(p,"tplot3b.png"))
    # tplot3c()
    # savefig(joinpath(p,"tplot3c.png"))
    tplot4()
    tplot5()
    splitf()
    splitf2()
    splitf3()
    demo()
    demo2()
    test_upper_env_dec()
end

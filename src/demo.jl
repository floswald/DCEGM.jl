function demos()
    m,p = runfdemo()
    d = Dict()
    d["T-1"] = plot(v_analytic(m,p,1,p.nT-1)[10:end], linecolor = :red, linewidth = 2, label = "work", legend = :bottomright, title = "Period T-2 Values",
    xlab = "cash on hand", ylab = "V(M)",ylims = (0,2.3))
    plot!(d["T-1"],v_analytic(m,p,2,p.nT-1)[103:200], linecolor = :blue, linewidth = 2, label = "retire")

    d["T-1x"] = plot(d["T-1"],xlims = (3,10),ylims = (1.25,1.5))
    d
end
"`plotPC(PCshock,horizon,G)`

Plots Philipps curve from MIT shock solutions for different shock sizes stored in `PCshock`
    `horizon` is the considered horizon in years (e.g. 1 year, 2 years,...)."
function plotPC(PCshock, horizon, G,legend)
    pPhilippsCurve = plot()
    pMCInflatCurve = plot()
    pOutptMCCurve = plot()
    prbinflCurve = plot()

    plotPC!(pPhilippsCurve, pMCInflatCurve, pOutptMCCurve, prbinflCurve, PCshock, horizon, G,legend,1,:circle)

    return pPhilippsCurve, pMCInflatCurve, pOutptMCCurve, prbinflCurve
end

"`plotPC!(p1,p2,p3,p4,PCshock,horizon,G)`

Adds Philipps curve plots from MIT shock solutions stored in `PCshock` to `p1`,`p2`,`p3`,`p4`
    `horizon` is the considered horizon in years (e.g. 1 year, 2 years,...)."
function plotPC!(p1, p2, p3, p4, PCshock, horizon, G, legend, color,shape)
    points = length(PCshock)
    PhilippsCurve = Array{Float64}(undef, points, 2)
    MCInflatCurve = Array{Float64}(undef, points, 2)
    OutptMCCurve = Array{Float64}(undef, points, 2)
    rbinflCurve = Array{Float64}(undef, points, 2)
    for i = 1:points
        IRF_horizon = findfirst(G.dt_cum .== minimum(G.dt_cum[G.dt_cum.>=horizon.*4.0]))

        PhilippsCurve[i, :] = [PCshock[i].Output_total_t[IRF_horizon]./PCshock[i].Output_total_t[end].-1.0, (PCshock[i].infl_t[IRF_horizon].-PCshock[i].infl_t[end]).*4.0].*100.0
        MCInflatCurve[i, :] = [PCshock[i].mc_t[IRF_horizon]./PCshock[i].mc_t[end].-1.0, (PCshock[i].infl_t[IRF_horizon].-PCshock[i].infl_t[end]).*4.0].*100.0
        OutptMCCurve[i, :] = [PCshock[i].Output_total_t[IRF_horizon]./PCshock[i].Output_total_t[end].-1.0, PCshock[i].mc_t[IRF_horizon]./PCshock[i].mc_t[end].-1.0].*100.0
        rbinflCurve[i, :] = [(PCshock[i].rb_t[IRF_horizon].-PCshock[i].rb_t[end]).*4.0, (PCshock[i].infl_t[IRF_horizon].-PCshock[i].infl_t[end]).*4.0].*100.0
    end
    plot!(p1, PhilippsCurve[:, 1], PhilippsCurve[:, 2], seriestype=:scatter, xlabel="Output gap (%)",ylabel="Inflation gap (pp. p.a.)",title="Phillips curve",legend =:topleft, label=legend, markershape = shape,markersize = 4, markerstrokewidth = 2,markerstrokecolor=color,color=color)
    plot!(p2, MCInflatCurve[:, 1], MCInflatCurve[:, 2], seriestype=:scatter, xlabel="Marginal cost gap (%)",ylabel="Inflation gap (pp. p.a.)",title="Marginal cost - inflation curve",legend =:topleft, label=legend, markershape = shape,markersize = 4, markerstrokewidth = 2,markerstrokecolor=color,color=color)
    plot!(p3, OutptMCCurve[:, 1], OutptMCCurve[:, 2], seriestype=:scatter, xlabel="Output gap (%)",ylabel="Marginal cost gap (%)",title="Output gap - marginal cost curve",legend =:topleft, label=legend, markershape = shape,markersize = 4, markerstrokewidth = 2,markerstrokecolor=color,color=color)
    plot!(p4, rbinflCurve[:, 1], rbinflCurve[:, 2], seriestype=:scatter, xlabel="Real interest rate gap (% p.a.)",ylabel="Inflation gap (pp. p.a.)",title="Real interest rate - inflation curve",legend =:topright, label=legend, markershape = shape,markersize = 4, markerstrokewidth = 2,markerstrokecolor=color,color=color)

end

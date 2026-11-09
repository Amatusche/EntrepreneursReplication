"`plotMITconv(MaxAbsCapDev,MaxAbsBondDev,dt_cum,K,Knew,it,ModelName)`

Plots convergence criteria over iterations in MIT shock solution and save as .png for inspection during computation."

function plotMITconv(MaxAbsCapDev,MaxAbsBondDev,dt_cum,K,Knew,it,ModelName; output_dir="../temp/")

    p1 = plot(1:it,MaxAbsCapDev[1:it],seriestype = :scatter,title="Max |Excess Capital Supply|",legend = false)
    p2 = plot(1:it,MaxAbsBondDev[1:it],seriestype = :scatter,title="Max |Excess Bond Supply|",legend = false)
    p3 = plot(dt_cum,Knew.-K,title="K_ new - K",legend=false)
    p = plot(p1,p2,p3,layout = (3,1), plot_title="Convergence")
    mkpath(output_dir)
    savefig(p, joinpath(output_dir, "IRFConvergence_"*string(ModelName)*".png"))

end

"`plotMIT!(plots,MIT,mp_shock_t,SSallend,years,G,ModelName)`

Plots IRFs over first `years` and save as png for inspection during computation."
function plotMIT!(plots,MIT,mp_shock_t,SSallend,years,G,ModelName; output_dir="../temp/")
    if G.dt_cum[end]>=years.*4.0
        IRF_horizon = findfirst(G.dt_cum.==minimum(G.dt_cum[G.dt_cum.>=years.*4.0]))
    else
        IRF_horizon = G.N_t
    end
    plots[1] = plot!(plots[1],G.dt_cum[1:IRF_horizon],4.0.*mp_shock_t[1:IRF_horizon])
    hline!(plots[1],0);
    title!(plots[1],"Shock (ann.)")

    plots[2] = plot!(plots[2],G.dt_cum[1:IRF_horizon],MIT.Output_total_t[1:IRF_horizon]./SSallend.Output_total.-1.0)
    hline!(plots[2],0);
    title!(plots[2],"Output")

    plots[3] = plot!(plots[3],G.dt_cum[1:IRF_horizon],MIT.AggrLabSup_t[1:IRF_horizon]./SSallend.AggrLabSup.-1.0)
    hline!(plots[3],0);;
    title!(plots[3],"Labor")

    plots[4] = plot!(plots[4],G.dt_cum[1:IRF_horizon],MIT.RealProfits_t[1:IRF_horizon]./SSallend.RealProfits.-1.0)
    hline!(plots[4],0);
    title!(plots[4],"Real profits")

    plots[5] = plot!(plots[5],G.dt_cum[1:IRF_horizon],MIT.Excess_Sup_Liq_t[1:IRF_horizon])
    hline!(plots[5],0);;
    title!(plots[5],"Excess liquid supply")

    plots[6] = plot!(plots[6],G.dt_cum[1:IRF_horizon],4.0.*MIT.rb_t[1:IRF_horizon].-4.0.*SSallend.rb)
    hline!(plots[6],0);;
    title!(plots[6],"rb (ann.)")

    plots[7] = plot!(plots[7],G.dt_cum[1:IRF_horizon],MIT.C_t[1:IRF_horizon]./SSallend.C.-1.0)
    hline!(plots[7],0);;
    title!(plots[7],"Consumption")

    plots[8] = plot!(plots[8],G.dt_cum[1:IRF_horizon],MIT.AggrCapSup_t[1:IRF_horizon]./SSallend.AggrCapSup.-1.0)
    hline!(plots[8],0);;
    title!(plots[8],"Market capital")

    plots[9] = plot!(plots[9],G.dt_cum[1:IRF_horizon],MIT.q_t[1:IRF_horizon]./SSallend.q.-1.0)
    hline!(plots[9],0);;
    title!(plots[9],"Stock price q (= market cap)")

    plots[11] = plot!(plots[11],G.dt_cum[1:IRF_horizon],4.0.*MIT.infl_t[1:IRF_horizon])
    hline!(plots[11],0);;
    title!(plots[11],"Infl. (ann.)")

    plots[12] = plot!(plots[12],G.dt_cum[1:IRF_horizon],MIT.Inv_t[1:IRF_horizon]./SSallend.Inv.-1.0)
    hline!(plots[12],0);
    title!(plots[12],"Investment")

    plots[13] = plot!(plots[13],G.dt_cum[1:IRF_horizon],MIT.InvC_M_t[1:IRF_horizon]./SSallend.InvC_M.-1.0)
    hline!(plots[13],0);
    title!(plots[13],"Investment (corp. market)")

    plots[14] = plot!(plots[14],G.dt_cum[1:IRF_horizon],MIT.Trans_t[1:IRF_horizon]./SSallend.Trans.-1.0)
    hline!(plots[14],0);;
    title!(plots[14],"Gov. transfer")

    plots[15] = plot!(plots[15],G.dt_cum[1:IRF_horizon],4.0.*MIT.ra_t[1:IRF_horizon].-4.0.*SSallend.ra)
    hline!(plots[15],0);;
    title!(plots[15],"ra (ann.)")

    plots[16] = plot!(plots[16],G.dt_cum[1:IRF_horizon],MIT.w_t[1:IRF_horizon]./SSallend.w.-1.0)
    hline!(plots[16],0);;
    title!(plots[16],"Real wage")

    plots[10] = plot!(plots[10],G.dt_cum[1:IRF_horizon],MIT.mc_t[1:IRF_horizon]./SSallend.pe.-1.0)
    hline!(plots[10],0);;
    title!(plots[10],"MC (= pe)")

    plots[17] = plot!(plots[17],G.dt_cum[1:IRF_horizon],MIT.InvE_t[1:IRF_horizon]./SSallend.InvE.-1.0)
    hline!(plots[17],0);
    title!(plots[17],"Investment (entre.)")

    plots[18] = plot!(plots[18],G.dt_cum[1:IRF_horizon],MIT.AggrFirmSav_t[1:IRF_horizon]./SSallend.AggrFirmSav.-1.0)
    hline!(plots[18],0);
    title!(plots[18],"Capital (private)")

    plots[19] = plot!(plots[19],G.dt_cum[1:IRF_horizon],MIT.rk_t[1:IRF_horizon].*4.0 .- SSallend.rk.*4.0)
    hline!(plots[19],0);
    title!(plots[19],"rk (ann.)")

    plots[20] = plot!(plots[20],G.dt_cum[1:IRF_horizon],MIT.tax_rev_entre_t[1:IRF_horizon]./SSallend.tax_rev_entre*100.0 .- 100.0)
    hline!(plots[20],0);
    title!(plots[20],"Tax revenue from entrepreneurs")

    plots[21] = plot!(plots[21],G.dt_cum[1:IRF_horizon],MIT.what_t[1:IRF_horizon]./SSallend.what.-1.0)
    hline!(plots[21],0);
    title!(plots[21],"Wage hat")

    plots[22] = plot!(plots[22],G.dt_cum[1:IRF_horizon],MIT.UnionProfits_t[1:IRF_horizon].-SSallend.UnionProfits)
    hline!(plots[22],0);
    title!(plots[22],"Union profits")

    plots[23] = plot!(plots[23],G.dt_cum[1:IRF_horizon],MIT.Tunion_t[1:IRF_horizon]./SSallend.Tunion.-1.0)
    hline!(plots[23],0);
    title!(plots[23],"Union lump sum taxes")

    plots[24] = plot!(plots[24],G.dt_cum[1:IRF_horizon],MIT.subs_payment_union_t[1:IRF_horizon]./SSallend.subs_payment_union.-1.0)
    hline!(plots[24],0);
    title!(plots[24],"Union subsidy")

    plots[25] = plot!(plots[25],G.dt_cum[1:IRF_horizon],4.0.*MIT.inflw_t[1:IRF_horizon])
    hline!(plots[25],0);;
    title!(plots[25],"Wage infl. (ann.)")

    plots[26] = plot!(plots[26], G.dt_cum[1:IRF_horizon-1],MIT.wagemarkdown_t[1:IRF_horizon-1]./SSallend.w)
    hline!(plots[26],0);
    title!(plots[26],"Wage markdown")

    plots[27] = plot(plots[27],G.dt_cum[1:IRF_horizon-1],MIT.dinflw_t[1:IRF_horizon-1]./ G.dt[1:IRF_horizon-1])
    hline!(plots[27],0);
    title!(plots[27],"d wage inflation")

    plots[28] = plot(plots[28],G.dt_cum[1:IRF_horizon-1],MIT.dLab_t[1:IRF_horizon-1]./ G.dt[1:IRF_horizon-1]./MIT.AggrLabSup_t[1:IRF_horizon-1])
    hline!(plots[28],0);
    title!(plots[28],"d N / N")

    plots[29] = plot(plots[29],G.dt_cum[1:IRF_horizon-1],diff(MIT.w_t[1:IRF_horizon])./ G.dt[1:IRF_horizon-1]./MIT.w_t[1:IRF_horizon-1])
    hline!(plots[29],0);
    title!(plots[29],"d w / w")

    p = plot(plots...,layout=(5,6),size=(1500,1500),legend=false)
    mkpath(output_dir)
    savefig(p, joinpath(output_dir, "IRFs"*string(years)*"years_"*string(ModelName)*".png"))

end

"`Statistics = calibstats!(SSall, Aentre, G,P,pHtm_MAT,wHtm_MAT,negLiq_MAT,pHtme_MAT,wHtme_MAT,negLiqe_MAT,simdebacker=false)`

Compute statistics for calibration form steady state and write to `Statistics`"
function calibstats!(SSall, Aentre, PolE, G, P, Opt, pHtm_MAT, wHtm_MAT, negLiq_MAT, pHtme_MAT, wHtme_MAT, negLiqe_MAT, simdebacker=false)

    pHtm_MAT .= ((G.b_grid_MAT .== 0.0) .| (G.b_grid_MAT .== P.underline_b)) .& (G.a_grid_MAT .== 0.0)
    wHtm_MAT .= ((G.b_grid_MAT .== 0.0) .| (G.b_grid_MAT .== P.underline_b)) .& (G.a_grid_MAT .> 0.0)
    negLiq_MAT .= G.b_grid_MAT .< 0.0
    if P.share_entre > 0
        pHtme_MAT .= ((G.b_grid_MAT_e .== 0.0) .| (G.b_grid_MAT_e .== P.underline_b)) .& (G.a_grid_MAT_e .== 0.0) .& (G.e_grid_MAT .== 0.0)
        wHtme_MAT .= ((G.b_grid_MAT_e .== 0.0) .| (G.b_grid_MAT_e .== P.underline_b)) .& ((G.a_grid_MAT_e .> 0.0) .| (G.e_grid_MAT .> 0.0))
        negLiqe_MAT .= G.b_grid_MAT_e .< 0.0
        pHtmShareW = pHtm_MAT[:]' * SSall.g[:]
        wHtmShareW = wHtm_MAT[:]' * SSall.g[:]
        pHtmShare = (pHtm_MAT[:]' * SSall.g[:]) * (1.0 .- P.share_entre) .+ (pHtme_MAT[:]' * SSall.ge[:]) .* P.share_entre
        wHtmShare = (wHtm_MAT[:]' * SSall.g[:]) * (1.0 .- P.share_entre) .+ (wHtme_MAT[:]' * SSall.ge[:]) .* P.share_entre
        htmShare = pHtmShare + wHtmShare
        negLiqShare = (negLiq_MAT[:]' * SSall.g[:]) .* (1.0 .- P.share_entre) .+ (negLiqe_MAT[:]' * SSall.ge[:]) .* P.share_entre
        LaborShareE = (SSall.AggrLabDemE .+ P.share_entre .* P.entre_lab) ./ (SSall.AggrLabSup .+ P.share_entre .* P.entre_lab)
        HtmShareE = pHtme_MAT[:]' * SSall.ge[:] .+ wHtme_MAT[:]' * SSall.ge[:]
        AggrIlliqSavAnnGDP = (SSall.AggrIlliqSav .+ SSall.AggrFirmSav) ./ (4.0 .* SSall.Output_total)
        AggrCapAnnGDP = (SSall.AggrIlliqSav .+ SSall.AggrFirmSav - SSall.q) ./ (4.0 .* SSall.Output_total)
        AggrLiqSavAnnGDP = SSall.AggrLiqSav ./ (4.0 * SSall.Output_total)
        ShareLiqE = SSall.AggrLiqSavE / SSall.AggrLiqSav
        ShareIlliqE = SSall.AggrIlliqSavE / SSall.AggrIlliqSav
        NetWorthMat = G.b_grid_MAT[:, :, 1, 1] .+ G.a_grid_MAT[:, :, 1, 1]
        NWsortIndex = sortperm(NetWorthMat[:])
        NW_grid = NetWorthMat[NWsortIndex]
        g_marg_NW_unsort = reshape(sum(SSall.g, dims=[4, 3]), G.N_b * G.N_a)
        g_marg_NW = g_marg_NW_unsort[NWsortIndex]
        Avg_NW_w = dot(g_marg_NW, NW_grid)
        # entrepreneurs
        NetWorthMat_e = G.b_grid_MAT_e[:, :, :, 1] .+ G.a_grid_MAT_e[:, :, :, 1] .+ G.e_grid_MAT[:, :, :, 1]
        NWsortIndex_e = sortperm(NetWorthMat_e[:])
        NW_grid_e = NetWorthMat_e[NWsortIndex_e]
        ge_marg_NW_unsort = reshape(sum(SSall.ge, dims=4), G.N_be * G.N_a * G.N_e)
        ge_marg_NW = ge_marg_NW_unsort[NWsortIndex_e]
        Avg_NW_e = dot(ge_marg_NW, NW_grid_e)
        Avg_NW_ratio = Avg_NW_e / Avg_NW_w
        Avg_Emp_Necessity = dot(SSall.LabDemE[:, :, :, 1][:], SSall.ge[:, :, :, 1][:]) ./ sum(SSall.ge[:, :, :, 1])
        Avg_Emp_Opportunity = dot(SSall.LabDemE[:, :, :, 2][:], SSall.ge[:, :, :, 2][:]) ./ sum(SSall.ge[:, :, :, 2])

        Emp_ratio = Avg_Emp_Necessity / Avg_Emp_Opportunity
        # DeBacker matrix: Form 1040 pass-through income, net of expected depreciation.
        if simdebacker

            rng = MersenneTwister(Opt.debacker_seed)
            _, DeBackerModelFull, _, _ = debackersim2(SSall.ge, Aentre, Opt.N_DeBackerSim1, SSall.TaxableBusinessIncomeE; rng=rng)
            mobility_stats = debacker_mobility_stats(DeBackerModelFull)
            Einc_immobility_ratio = mobility_stats.immobility_ratio
            EInc_rank_autocorrelation = mobility_stats.rank_autocorrelation
            EInc_large_loss_probability = mobility_stats.large_loss_probability
            EInc_large_gain_probability = mobility_stats.large_gain_probability
            EInc_avg_deciles_lost = mobility_stats.avg_deciles_lost
        else
            Einc_immobility_ratio = NaN
            EInc_rank_autocorrelation = NaN
            EInc_large_loss_probability = NaN
            EInc_large_gain_probability = NaN
            EInc_avg_deciles_lost = NaN
        end

    elseif P.share_entre == 0
        pHtmShare = pHtm_MAT[:]' * SSall.g[:]
        wHtmShare = wHtm_MAT[:]' * SSall.g[:]
        pHtmShareW = pHtmShare
        wHtmShareW = wHtmShare
        htmShare = pHtmShare .+ wHtmShare
        negLiqShare = negLiq_MAT[:]' * SSall.g[:]
        AggrIlliqSavAnnGDP = SSall.AggrIlliqSav ./ (4.0 .* SSall.Output_total)
        AggrLiqSavAnnGDP = SSall.AggrLiqSav ./ (4.0 .* SSall.Output_total)
        Avg_NW_ratio = 0.0
        Emp_ratio = NaN
    end

    Statistics = CalibTargets(AggrIlliqSavAnnGDP, AggrLiqSavAnnGDP,
        pHtmShare, wHtmShare, negLiqShare, LaborShareE, ShareLiqE, ShareIlliqE, HtmShareE, Avg_NW_ratio, AggrCapAnnGDP, Emp_ratio, Einc_immobility_ratio, EInc_rank_autocorrelation, EInc_large_loss_probability, EInc_large_gain_probability, EInc_avg_deciles_lost)

end

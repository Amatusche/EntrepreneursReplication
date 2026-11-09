" `stats(SSall, PolE, Aentre, G, P, Opt)`

Compute various statistics of interest for the steady state.

"

function stats(SSall, PolE, Aentre, G, P, Opt)
    println("Computing steady state statistics")

    # Targeted moments
    pHtm_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    wHtm_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    negLiq_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    pHtme_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    wHtme_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    negLiqe_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    CalibStats = calibstats!(SSall, Aentre, PolE, G, P, Opt, pHtm_MAT, wHtm_MAT, negLiq_MAT, pHtme_MAT, wHtme_MAT, negLiqe_MAT, false)
    @unpack AggrIlliqSavAnnGDP, AggrLiqSavAnnGDP,
    pHtmShare, wHtmShare, negLiqShare, LaborShareE, ShareLiqE, ShareIlliqE, HtmShareE, Emp_ratio = CalibStats

    θ = reshape(PolE.θ, G.N_be, G.N_a, G.N_e, G.N_y)

    # worker
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
    NWshareE = Avg_NW_e * P.share_entre / (Avg_NW_e * P.share_entre + Avg_NW_w * (1.0 - P.share_entre))
    # joint
    NWsortIndex_joint = sortperm([NW_grid_e; NW_grid])
    NW_grid_joint = [NW_grid_e; NW_grid][NWsortIndex_joint]
    g_marg_NW_joint = [P.share_entre * ge_marg_NW; (1 - P.share_entre) * g_marg_NW]
    g_marg_NW_joint = g_marg_NW_joint[NWsortIndex_joint]

    # Untargeted moments I: Wealth distribution
    NWShares = Array{Float64}(undef, 7, 2)
    NWShares[:, 1] .= [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001]
    NWShares[:, 2] .= [percshare(NW_grid_joint, g_marg_NW_joint, p) for p in NWShares[:, 1]]

    # Untargetd moments II: Entrepreneurs by net worth quantiles
    ShareEByNWGroup = Array{Float64}(undef, 19, 3)
    ShareEByNWGroup[:, 1] .= [0.0, 0.0, 0.4, 0.6, 0.8, 0.9, 0.95, 0.98, 0.99, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # lower bounds
    ShareEByNWGroup[:, 2] .= [1.0, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # upper bounds
    NW_perc_up = Array{Float64}(undef, length(ShareEByNWGroup[:, 1]))
    NW_perc_lo = Array{Float64}(undef, length(ShareEByNWGroup[:, 1]))
    ge_partconsidered = repeat(SSall.ge_marg_bae[:], 1, length(ShareEByNWGroup[:, 2]))
    g_partconsidered = repeat(SSall.g_marg_ba[:], 1, length(ShareEByNWGroup[:, 2]))
    for gg = eachindex(ShareEByNWGroup[:, 1])
        NW_perc_up[gg] = NW_grid_joint[findlast(cumsum(g_marg_NW_joint) .< ShareEByNWGroup[gg, 2])]
        NW_perc_lo[gg] = NW_grid_joint[findfirst(cumsum(g_marg_NW_joint) .>= ShareEByNWGroup[gg, 1])]
        # crop distributions to part between respective cutoffs
        g_partconsidered[((NetWorthMat[:].>=NW_perc_lo[gg]).&(NetWorthMat[:].<NW_perc_up[gg])).==0, gg] .= 0
        ge_partconsidered[((NetWorthMat_e[:].>=NW_perc_lo[gg]).&(NetWorthMat_e[:].<NW_perc_up[gg])).==0, gg] .= 0
        ShareEByNWGroup[gg, 3] = sum(ge_partconsidered[:, gg]) .* P.share_entre ./ (sum(ge_partconsidered[:, gg]) .* P.share_entre .+ sum(g_partconsidered[:, gg]) .* (1.0 .- P.share_entre))
    end

    # Simulate annual taxable income for the DeBacker mobility matrix and compute
    # business returns as annual taxable business income divided by business equity.

    # Cash-flow returns subtract net investment and adjustment costs from taxable income.
    @warn "DeBack Nsim = " * string(Opt.N_DeBackerSim1)
    rng = MersenneTwister(Opt.debacker_seed)
    DeBackerModelQuint, DeBackerModelFull, Income, states_current = debackersim2(SSall.ge, Aentre, Opt.N_DeBackerSim1, SSall.TaxableBusinessIncomeE; rng=rng)
    immobility_ratio = mean(diag(DeBackerModelQuint))

    firm_adjustment_cost = P.χe1 .* abs.(PolE.f ./ G.e_grid_MAT_adj[:]) .^ P.χe2 .* G.e_grid_MAT_adj[:]
    CashFlowBusinessIncomeE = SSall.TaxableBusinessIncomeE .-
                              reshape(PolE.f .+ firm_adjustment_cost, size(SSall.TaxableBusinessIncomeE))
    rng_cashflow = MersenneTwister(Opt.debacker_seed)
    _, _, IncomeCF, states_current_cf = debackersim2(SSall.ge, Aentre, Opt.N_DeBackerSim1, CashFlowBusinessIncomeE; rng=rng_cashflow)
    @assert states_current_cf == states_current "Return simulations with the same seed produced different state paths"

    MedReturn = Array{Float64}(undef, 10, 3)
    MedReturn[:, 1] .= collect(0:0.1:0.9)
    MedReturn[:, 2] .= MedReturn[:, 1] .+ 0.1
    MedReturnCF = similar(MedReturn)
    MedReturnCF[:, 1:2] .= MedReturn[:, 1:2]
    NetWorthMat_e_full = G.b_grid_MAT_e[states_current] .+ G.a_grid_MAT_e[states_current] .+ G.e_grid_MAT[states_current]
    business_equity = G.e_grid_MAT[states_current]
    busreturn = Income ./ business_equity
    busreturnCF = IncomeCF ./ business_equity
    for gg = eachindex(MedReturn[:, 1])
        ge_partconsidered = (NetWorthMat_e_full .>= quantile(NetWorthMat_e_full, MedReturn[gg, 1])) .& (NetWorthMat_e_full .< quantile(NetWorthMat_e_full, MedReturn[gg, 2])) .& (business_equity .> 0.0)
        MedReturn[gg, 3] = any(ge_partconsidered) ? median(busreturn[ge_partconsidered]) : NaN
        MedReturnCF[gg, 3] = any(ge_partconsidered) ? median(busreturnCF[ge_partconsidered]) : NaN
    end

    # TFP
    logTFPE = log(SSall.AggrOutputE) .-
              (1.0 - P.alphaaE) .* P.nu .* log(SSall.AggrLabDemE) .-
              P.alphaaE .* P.nu .* log(SSall.AggrFirmSav + SSall.AggrCapDemE)
    logTFPC = log(SSall.OutputC) .- (1.0 - P.alphaa) .* log(SSall.AggrLabDemC) .- P.alphaa .* log(SSall.AggrCapDemC)

    #  Share of capital at private business
    PrivateFirmCapShare = (SSall.AggrFirmSav + SSall.AggrCapDemE) / (SSall.AggrFirmSav + SSall.AggrCapDemE + SSall.AggrCapDemC)

    # Firmsize distribution
    # Sales
    Sales = Array{Float64}(undef, 4, 4)
    Sales[:, 1] .= [0.0, 2.0, 10.0, 50.0] .* 1000000.0
    Sales[:, 2] .= [2.0, 10.0, 50.0, Inf] .* 1000000.0
    SalesMat = SSall.pe .* SSall.OutputE .* P.NomGDPFac * 4.0
    Index = sortperm(SalesMat[:])
    SalesMat = SalesMat[Index]
    FirmSizeDistrNo = SSall.ge[Index]
    FirmSizeDistrSales = FirmSizeDistrNo .* SalesMat
    FirmSizeDistrSales = FirmSizeDistrSales ./ sum(FirmSizeDistrSales)
    Sales[:, 3] = [sum(FirmSizeDistrNo[(SalesMat.>=Sales[gg, 1]).&(SalesMat.<Sales[gg, 2])]) for gg in eachindex(Sales[:, 1])]
    Sales[:, 4] = [sum(FirmSizeDistrSales[(SalesMat.>=Sales[gg, 1]).&(SalesMat.<Sales[gg, 2])]) for gg in eachindex(Sales[:, 1])]

    # Employment
    Employment = Array{Float64}(undef, 5, 4)
    Employment[:, 1] .= [0.0, 1.0, 10.0, 50, 250.0]
    Employment[:, 2] .= [1.0, 10.0, 50.0, 250.0, Inf]
    EmploymentMat = SSall.LabDemE[:]
    EmploymentMat = EmploymentMat ./ SSall.AggrLabSup
    Index = sortperm(EmploymentMat[:])
    EmploymentMat = EmploymentMat[Index]
    FirmSizeDistrNo = SSall.ge[Index]

    FirmSizeDistrEmployment = FirmSizeDistrNo .* EmploymentMat
    FirmSizeDistrEmployment = FirmSizeDistrEmployment ./ sum(FirmSizeDistrEmployment)
    Employment[:, 3] = [sum(FirmSizeDistrNo[(EmploymentMat.>=Employment[gg, 1]).&(EmploymentMat.<Employment[gg, 2])]) for gg in eachindex(Employment[:, 1])]
    Employment[:, 4] = [sum(FirmSizeDistrEmployment[(EmploymentMat.>=Employment[gg, 1]).&(EmploymentMat.<Employment[gg, 2])]) for gg in eachindex(Employment[:, 1])]

    Stats = SteadyStats(AggrIlliqSavAnnGDP, AggrLiqSavAnnGDP, pHtmShare, wHtmShare, negLiqShare, LaborShareE, PrivateFirmCapShare, ShareLiqE, ShareIlliqE, HtmShareE, MedReturn, MedReturnCF, ShareEByNWGroup, NWShares, DeBackerModelQuint, DeBackerModelFull, Avg_NW_ratio, Emp_ratio)

    return Stats::SteadyStats

end

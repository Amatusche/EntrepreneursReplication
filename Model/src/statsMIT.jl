" `statsMIT(MIT, TransPolW, TransPolE, MITgenerators, Partial, Partialgenerators, SSall, Pol, PolE, A, Aentre, G, P, Opt)

Compute various statistics of interest for the MIT shock.

"

function statsMIT(MIT, TransPolW, TransPolE, MITgenerators, Partial, Partialgenerators, SSall, Pol, PolE, A, Aentre, G, P, Opt=options())
    println("Computing MIT statistics")

    Δrb = Partial.MIT_partial.rb_t[1] - SSall.rb
    AggrPartialSemiElaInvE_nonM = (Partial.MIT_partial.InvE_nonM_t[1] ./ SSall.InvE_nonM .- 1.0) ./ Δrb
    AggrPartialSemiElaInv_M = (Partial.MIT_partial.Inv_M_t[1] ./ SSall.Inv_M .- 1.0) ./ Δrb
    AggrPartialSemiElaInvCW = (Partial.MIT_partial.InvCW_t[1] ./ SSall.InvCW .- 1.0) ./ Δrb
    AggrPartialSemiElaInvCE = (Partial.MIT_partial.InvCE_t[1] ./ SSall.InvCE .- 1.0) ./ Δrb
    AggrPartialSemiElaInvE = (Partial.MIT_partial.InvE_t[1] ./ SSall.InvE .- 1.0) ./ Δrb
    AggrPartialSemiElaInvC_M = (Partial.MIT_partial.InvC_M_t[1] ./ SSall.InvC_M .- 1.0) ./ Δrb
    AggrPartialSemiElaInv = (Partial.MIT_partial.Inv_t[1] ./ SSall.Inv .- 1.0) ./ Δrb

    # Compute htm, negative liquidy and leverage matrices
    pHtm_MAT = ((G.b_grid_MAT .== 0.0) .| (G.b_grid_MAT .== P.underline_b)) .& (G.a_grid_MAT .== 0.0)
    wHtm_MAT = ((G.b_grid_MAT .== 0.0) .| (G.b_grid_MAT .== P.underline_b)) .& (G.a_grid_MAT .> 0.0)
    negLiq_MAT = G.b_grid_MAT .< 0.0
    pHtm_MAT_e = ((G.b_grid_MAT_e .== 0.0) .| (G.b_grid_MAT_e .== P.underline_be)) .& (G.a_grid_MAT_e .== 0.0) .& (G.e_grid_MAT .== 0.0)
    wHtm_MAT_e = ((G.b_grid_MAT_e .== 0.0) .| (G.b_grid_MAT_e .== P.underline_be)) .& ((G.a_grid_MAT_e .> 0.0) .| (G.e_grid_MAT .> 0.0))
    negLiq_MAT_e = G.b_grid_MAT_e .< 0.0
    θ_MAT = reshape(PolE.θ, G.N_be, G.N_a, G.N_e, G.N_y)
    PrivateFirmCap_MAT = θ_MAT .* G.e_grid_MAT

    NetWorthMat_e_full = G.b_grid_MAT_e + G.a_grid_MAT_e + G.e_grid_MAT
    NWsortIndex_e_full = sortperm(NetWorthMat_e_full[:])
    NW_grid_e_full = NetWorthMat_e_full[NWsortIndex_e_full]
    ge_NW_unsort = SSall.ge[:]
    ge_NW_sort = ge_NW_unsort[NWsortIndex_e_full]

    NetWorthMat_e = NetWorthMat_e_full[:, :, :, 1]
    NWsortIndex_e = sortperm(vec(NetWorthMat_e))
    NW_grid_e = NetWorthMat_e[NWsortIndex_e]
    NWgridsmooth = range(minimum(NW_grid_e), stop=maximum(NW_grid_e), length=1000)
    ge1_NW_unsort = SSall.ge[:, :, :, 1][:] ./ sum(SSall.ge[:, :, :, 1])
    ge1_NW_sort = ge1_NW_unsort[NWsortIndex_e]
    ge2_NW_unsort = SSall.ge[:, :, :, 2][:] ./ sum(SSall.ge[:, :, :, 2])
    ge2_NW_sort = ge2_NW_unsort[NWsortIndex_e]

    # Estimate elasticities and MPIs over net worth distribution using local constant regressions
    steps1 = 0.1
    steps2 = 0.02
    startfine = 0.9
    bounds_lo = vcat(0:steps1:startfine, (startfine+steps2):steps2:(1-steps2))
    bounds_up = vcat(bounds_lo[2:end], 1)

    d_ela_MAT_MIT = reshape((TransPolW[1].d .- Pol.d), G.N_b, G.N_a, G.N_j, G.N_k) ./ G.a_grid_MAT_adj
    de_ela_MAT_MIT = reshape((TransPolE[1].de .- PolE.de), G.N_be, G.N_a, G.N_e, G.N_y) ./ G.a_grid_MAT_e_adj
    f_ela_MAT_MIT = reshape((TransPolE[1].f .- PolE.f), G.N_be, G.N_a, G.N_e, G.N_y) ./ G.e_grid_MAT_adj
    θ_ela_MAT_MIT = reshape((TransPolE[1].θ .- PolE.θ), G.N_be, G.N_a, G.N_e, G.N_y)

    f_ela_avg = dot(SSall.ge[:], f_ela_MAT_MIT[:]) / sum(SSall.ge[:])
    θ_ela_avg = dot(SSall.ge[:], θ_ela_MAT_MIT[:]) / sum(SSall.ge[:])
    d_ela_avg = dot(SSall.g[:], d_ela_MAT_MIT[:]) / sum(SSall.g[:])
    de_ela_avg = dot(SSall.ge[:], de_ela_MAT_MIT[:]) / sum(SSall.ge[:])
    d_ela_total_avg = (1.0 - P.share_entre) * d_ela_avg + P.share_entre * de_ela_avg

    f_partial_ela_avg = dot(SSall.ge[:], Partial.f_ela_MAT[:]) / sum(SSall.ge[:])
    θ_partial_ela_avg = dot(SSall.ge[:], Partial.θ_ela_MAT[:]) / sum(SSall.ge[:])
    d_partial_ela_avg = dot(SSall.g[:], Partial.d_ela_MAT[:]) / sum(SSall.g[:])
    de_partial_ela_avg = dot(SSall.ge[:], Partial.de_ela_MAT[:]) / sum(SSall.ge[:])
    d_partial_ela_total_avg = (1.0 - P.share_entre) * d_partial_ela_avg + P.share_entre * de_partial_ela_avg

    # Low productivity entrepreneurs
    ge1_partconsidered, perc_up1, perc_lo1 = percmask(SSall.ge[:, :, :, 1], bounds_lo, bounds_up, NetWorthMat_e)

    ElaMIT_f_OverNWPerc1 = distributionavg(f_ela_MAT_MIT[:, :, :, 1], ge1_partconsidered)
    ElaMIT_de_OverNWPerc1 = distributionavg(de_ela_MAT_MIT[:, :, :, 1], ge1_partconsidered)
    ElaMIT_θ_OverNWPerc1 = distributionavg(θ_ela_MAT_MIT[:, :, :, 1], ge1_partconsidered)
    ElaPartial_f_OverNWPerc1 = distributionavg(Partial.f_ela_MAT[:, :, :, 1], ge1_partconsidered)
    ElaPartial_de_OverNWPerc1 = distributionavg(Partial.de_ela_MAT[:, :, :, 1], ge1_partconsidered)
    ElaPartial_θ_OverNWPerc1 = distributionavg(Partial.θ_ela_MAT[:, :, :, 1], ge1_partconsidered)

    PrivateFirmOwnCap1 = distributionavg(G.e_grid_MAT[:, :, :, 1], ge1_partconsidered) .* vec(sum(ge1_partconsidered, dims=1)) .* P.share_entre ./ (perc_up1 .- perc_lo1)
    PrivateFirmCap1 = distributionavg(PrivateFirmCap_MAT[:, :, :, 1], ge1_partconsidered) .* vec(sum(ge1_partconsidered, dims=1)) .* P.share_entre ./ (perc_up1 .- perc_lo1)
    IlliqAss_e1 = distributionavg(G.a_grid_MAT_e[:, :, :, 1], ge1_partconsidered) .* vec(sum(ge1_partconsidered, dims=1)) .* P.share_entre ./ (perc_up1 .- perc_lo1)
    LiqAss_e1 = distributionavg(G.b_grid_MAT_e[:, :, :, 1], ge1_partconsidered) .* vec(sum(ge1_partconsidered, dims=1)) .* P.share_entre ./ (perc_up1 .- perc_lo1)
    wHtmShare_e1 = distributionavg(wHtm_MAT_e[:, :, :, 1], ge1_partconsidered)
    pHtmShare_e1 = distributionavg(pHtm_MAT_e[:, :, :, 1], ge1_partconsidered)
    negLiqShare_e1 = distributionavg(negLiq_MAT_e[:, :, :, 1], ge1_partconsidered)

    HtmShare_e1 = wHtmShare_e1 + pHtmShare_e1

    nodes1 = 0.5 * (perc_up1 + perc_lo1)
    itp = CubicSplineInterpolator(nodes1, ElaMIT_f_OverNWPerc1, NoBoundaries())
    ElaMIT_f_smooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, ElaMIT_θ_OverNWPerc1, NoBoundaries())
    ElaMIT_θ_smooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, ElaMIT_de_OverNWPerc1, NoBoundaries())
    ElaMIT_de_smooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, ElaPartial_f_OverNWPerc1, NoBoundaries())
    ElaPartial_f_smooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, ElaPartial_θ_OverNWPerc1, NoBoundaries())
    ElaPartial_θ_smooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, ElaPartial_de_OverNWPerc1, NoBoundaries())
    ElaPartial_de_smooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, PrivateFirmOwnCap1, NoBoundaries())
    PrivateFirmOwnCapSmooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, PrivateFirmCap1, NoBoundaries())
    PrivateFirmCapSmooth1 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes1, IlliqAss_e1, NoBoundaries())
    IlliqAss_e_smooth1 = itp.(NWgridsmooth)

    # High productivity entrepreneurs
    ge2_partconsidered, perc_up2, perc_lo2 = percmask(SSall.ge[:, :, :, 2], bounds_lo, bounds_up, NetWorthMat_e)

    ElaMIT_f_OverNWPerc2 = distributionavg(f_ela_MAT_MIT[:, :, :, 2], ge2_partconsidered)
    ElaMIT_de_OverNWPerc2 = distributionavg(de_ela_MAT_MIT[:, :, :, 2], ge2_partconsidered)
    ElaMIT_θ_OverNWPerc2 = distributionavg(θ_ela_MAT_MIT[:, :, :, 2], ge2_partconsidered)
    ElaPartial_f_OverNWPerc2 = distributionavg(Partial.f_ela_MAT[:, :, :, 2], ge2_partconsidered)
    ElaPartial_de_OverNWPerc2 = distributionavg(Partial.de_ela_MAT[:, :, :, 2], ge2_partconsidered)
    ElaPartial_θ_OverNWPerc2 = distributionavg(Partial.θ_ela_MAT[:, :, :, 2], ge2_partconsidered)

    PrivateFirmOwnCap2 = distributionavg(G.e_grid_MAT[:, :, :, 2], ge2_partconsidered) .* vec(sum(ge2_partconsidered, dims=1)) .* P.share_entre ./ (perc_up2 .- perc_lo2)
    PrivateFirmCap2 = distributionavg(PrivateFirmCap_MAT[:, :, :, 2], ge2_partconsidered) .* vec(sum(ge2_partconsidered, dims=1)) .* P.share_entre ./ (perc_up2 .- perc_lo2)
    IlliqAss_e2 = distributionavg(G.a_grid_MAT_e[:, :, :, 2], ge2_partconsidered) .* vec(sum(ge2_partconsidered, dims=1)) .* P.share_entre ./ (perc_up2 .- perc_lo2)
    LiqAss_e2 = distributionavg(G.b_grid_MAT_e[:, :, :, 2], ge2_partconsidered) .* vec(sum(ge2_partconsidered, dims=1)) .* P.share_entre ./ (perc_up2 .- perc_lo2)
    wHtmShare_e2 = distributionavg(wHtm_MAT_e[:, :, :, 2], ge2_partconsidered)
    pHtmShare_e2 = distributionavg(pHtm_MAT_e[:, :, :, 2], ge2_partconsidered)
    negLiqShare_e2 = distributionavg(negLiq_MAT_e[:, :, :, 2], ge2_partconsidered)

    HtmShare_e2 = wHtmShare_e2 + pHtmShare_e2

    nodes2 = 0.5 * (perc_up2 + perc_lo2)
    itp = CubicSplineInterpolator(nodes2, ElaMIT_f_OverNWPerc2, NoBoundaries())
    ElaMIT_f_smooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, ElaMIT_θ_OverNWPerc2, NoBoundaries())
    ElaMIT_θ_smooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, ElaMIT_de_OverNWPerc2, NoBoundaries())
    ElaMIT_de_smooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, ElaPartial_f_OverNWPerc2, NoBoundaries())
    ElaPartial_f_smooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, ElaPartial_θ_OverNWPerc2, NoBoundaries())
    ElaPartial_θ_smooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, ElaPartial_de_OverNWPerc2, NoBoundaries())
    ElaPartial_de_smooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, PrivateFirmOwnCap2, NoBoundaries())
    PrivateFirmOwnCapSmooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, PrivateFirmCap2, NoBoundaries())
    PrivateFirmCapSmooth2 = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes2, IlliqAss_e2, NoBoundaries())
    IlliqAss_e_smooth2 = itp.(NWgridsmooth)

    # All entrepreneurs jointly
    ge_partconsidered, perc_up, perc_lo = percmask(SSall.ge, bounds_lo, bounds_up, NetWorthMat_e_full)

    ElaMIT_f_OverNWPerc = distributionavg(f_ela_MAT_MIT, ge_partconsidered)
    ElaMIT_de_OverNWPerc = distributionavg(de_ela_MAT_MIT, ge_partconsidered)
    ElaMIT_θ_OverNWPerc = distributionavg(θ_ela_MAT_MIT, ge_partconsidered)
    ElaPartial_f_OverNWPerc = distributionavg(Partial.f_ela_MAT, ge_partconsidered)
    ElaPartial_de_OverNWPerc = distributionavg(Partial.de_ela_MAT, ge_partconsidered)
    ElaPartial_θ_OverNWPerc = distributionavg(Partial.θ_ela_MAT, ge_partconsidered)

    PrivateFirmOwnCap = distributionavg(G.e_grid_MAT, ge_partconsidered) .* vec(sum(ge_partconsidered, dims=1)) .* P.share_entre ./ (perc_up .- perc_lo)
    PrivateFirmCap = distributionavg(PrivateFirmCap_MAT, ge_partconsidered) .* vec(sum(ge_partconsidered, dims=1)) .* P.share_entre ./ (perc_up .- perc_lo)
    IlliqAss_e = distributionavg(G.a_grid_MAT_e, ge_partconsidered) .* vec(sum(ge_partconsidered, dims=1)) .* P.share_entre ./ (perc_up .- perc_lo)
    LiqAss_e = distributionavg(G.b_grid_MAT_e, ge_partconsidered) .* vec(sum(ge_partconsidered, dims=1)) .* P.share_entre ./ (perc_up .- perc_lo)
    wHtmShare_e = distributionavg(wHtm_MAT_e, ge_partconsidered)
    pHtmShare_e = distributionavg(pHtm_MAT_e, ge_partconsidered)
    negLiqShare_e = distributionavg(negLiq_MAT_e, ge_partconsidered)

    HtmShare_e = wHtmShare_e + pHtmShare_e

    NWgridsmooth = range(minimum(NW_grid_e), stop=maximum(NW_grid_e), length=1000)
    nodes = 0.5 * (perc_up + perc_lo)
    itp = CubicSplineInterpolator(nodes, ElaMIT_f_OverNWPerc, NoBoundaries())
    ElaMIT_f_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, ElaMIT_θ_OverNWPerc, NoBoundaries())
    ElaMIT_θ_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, ElaMIT_de_OverNWPerc, NoBoundaries())
    ElaMIT_de_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, ElaPartial_f_OverNWPerc, NoBoundaries())
    ElaPartial_f_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, ElaPartial_θ_OverNWPerc, NoBoundaries())
    ElaPartial_θ_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, ElaPartial_de_OverNWPerc, NoBoundaries())
    ElaPartial_de_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, PrivateFirmOwnCap, NoBoundaries())
    PrivateFirmOwnCapSmooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, PrivateFirmCap, NoBoundaries())
    PrivateFirmCapSmooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodes, IlliqAss_e, NoBoundaries())
    IlliqAss_e_smooth = itp.(NWgridsmooth)

    # Workers
    NetWorthMat_full = G.b_grid_MAT + G.a_grid_MAT
    g_partconsidered, perc_upw, perc_low = percmask(SSall.g, bounds_lo, bounds_up, NetWorthMat_full)

    ElaMIT_d_OverNWPerc = distributionavg(d_ela_MAT_MIT, g_partconsidered)
    ElaPartial_d_OverNWPerc = distributionavg(Partial.d_ela_MAT, g_partconsidered)

    IlliqAss = distributionavg(G.a_grid_MAT, g_partconsidered) .* vec(sum(g_partconsidered, dims=1)) .* (1.0 .- P.share_entre) ./ (perc_upw .- perc_low)
    LiqAss = distributionavg(G.b_grid_MAT, g_partconsidered) .* vec(sum(g_partconsidered, dims=1)) .* (1.0 .- P.share_entre) ./ (perc_upw .- perc_low)
    wHtmShare = distributionavg(wHtm_MAT, g_partconsidered)
    pHtmShare = distributionavg(pHtm_MAT, g_partconsidered)
    negLiqShare = distributionavg(negLiq_MAT, g_partconsidered)

    HtmShare = wHtmShare + pHtmShare

    nodesw = 0.5 * (perc_upw + perc_low)
    itp = CubicSplineInterpolator(nodesw, ElaMIT_d_OverNWPerc, NoBoundaries())
    ElaMIT_d_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodesw, ElaPartial_d_OverNWPerc, NoBoundaries())
    ElaPartial_d_smooth = itp.(NWgridsmooth)
    itp = CubicSplineInterpolator(nodesw, IlliqAss, NoBoundaries())
    IlliqAss_smooth = itp.(NWgridsmooth)

    # Common wealth groups used in the bequest capital-distribution figure.
    nw_capital_shares = capital_shares_by_net_worth(SSall, PolE, G, P)
    PrivateFirmOwnCapByNWShare = nw_capital_shares.PrivateFirmOwnCapByNWShare
    PrivateFirmCapByNWShare = nw_capital_shares.PrivateFirmCapByNWShare
    IlliqAssShare = nw_capital_shares.IlliqAssShare

    MITStats = (NWgridsmooth=NWgridsmooth,
        nodes1=nodes1,
        ElaMIT_f_OverNWPerc1=ElaMIT_f_OverNWPerc1,
        ElaMIT_de_OverNWPerc1=ElaMIT_de_OverNWPerc1,
        ElaPartial_f_OverNWPerc1=ElaPartial_f_OverNWPerc1,
        ElaPartial_de_OverNWPerc1=ElaPartial_de_OverNWPerc1,
        PrivateFirmOwnCap1=PrivateFirmOwnCap1,
        PrivateFirmCap1=PrivateFirmCap1,
        IlliqAss_e1=IlliqAss_e1,
        LiqAss_e1=LiqAss_e1,
        wHtmShare_e1=wHtmShare_e1,
        pHtmShare_e1=pHtmShare_e1,
        negLiqShare_e1=negLiqShare_e1,
        HtmShare_e1=HtmShare_e1,
        ElaMIT_f_smooth1=ElaMIT_f_smooth1,
        ElaMIT_de_smooth1=ElaMIT_de_smooth1,
        ElaPartial_f_smooth1=ElaPartial_f_smooth1,
        ElaPartial_de_smooth1=ElaPartial_de_smooth1,
        PrivateFirmOwnCapSmooth1=PrivateFirmOwnCapSmooth1,
        PrivateFirmCapSmooth1=PrivateFirmCapSmooth1,
        IlliqAss_e_smooth1=IlliqAss_e_smooth1, nodes2=nodes2,
        ElaMIT_f_OverNWPerc2=ElaMIT_f_OverNWPerc2,
        ElaMIT_de_OverNWPerc2=ElaMIT_de_OverNWPerc2,
        ElaPartial_f_OverNWPerc2=ElaPartial_f_OverNWPerc2,
        ElaPartial_de_OverNWPerc2=ElaPartial_de_OverNWPerc2,
        PrivateFirmOwnCap2=PrivateFirmOwnCap2,
        PrivateFirmCap2=PrivateFirmCap2,
        IlliqAss_e2=IlliqAss_e2,
        LiqAss_e2=LiqAss_e2,
        wHtmShare_e2=wHtmShare_e2,
        pHtmShare_e2=pHtmShare_e2,
        negLiqShare_e2=negLiqShare_e2,
        HtmShare_e2=HtmShare_e2,
        ElaMIT_f_smooth2=ElaMIT_f_smooth2,
        ElaMIT_de_smooth2=ElaMIT_de_smooth2,
        ElaPartial_f_smooth2=ElaPartial_f_smooth2,
        ElaPartial_de_smooth2=ElaPartial_de_smooth2,
        PrivateFirmOwnCapSmooth2=PrivateFirmOwnCapSmooth2,
        PrivateFirmCapSmooth2=PrivateFirmCapSmooth2,
        IlliqAss_e_smooth2=IlliqAss_e_smooth2, nodes=nodes,
        ElaMIT_f_OverNWPerc=ElaMIT_f_OverNWPerc,
        ElaMIT_de_OverNWPerc=ElaMIT_de_OverNWPerc,
        ElaPartial_f_OverNWPerc=ElaPartial_f_OverNWPerc,
        ElaPartial_de_OverNWPerc=ElaPartial_de_OverNWPerc,
        PrivateFirmOwnCap=PrivateFirmOwnCap,
        PrivateFirmCap=PrivateFirmCap,
        IlliqAss_e=IlliqAss_e,
        LiqAss_e=LiqAss_e,
        wHtmShare_e=wHtmShare_e,
        pHtmShare_e=pHtmShare_e,
        negLiqShare_e=negLiqShare_e,
        HtmShare_e=HtmShare_e,
        ElaMIT_f_smooth=ElaMIT_f_smooth,
        ElaMIT_de_smooth=ElaMIT_de_smooth,
        ElaPartial_f_smooth=ElaPartial_f_smooth,
        ElaPartial_de_smooth=ElaPartial_de_smooth,
        PrivateFirmOwnCapSmooth=PrivateFirmOwnCapSmooth,
        PrivateFirmCapSmooth=PrivateFirmCapSmooth,
        IlliqAss_e_smooth=IlliqAss_e_smooth, nodesw=nodesw, ElaMIT_d_OverNWPerc=ElaMIT_d_OverNWPerc,
        ElaPartial_d_OverNWPerc=ElaPartial_d_OverNWPerc,
        IlliqAss=IlliqAss,
        LiqAss=LiqAss,
        wHtmShare=wHtmShare,
        pHtmShare=pHtmShare,
        negLiqShare=negLiqShare,
        HtmShare=HtmShare,
        ElaMIT_d_smooth=ElaMIT_d_smooth,
        ElaPartial_d_smooth=ElaPartial_d_smooth,
        IlliqAss_smooth=IlliqAss_smooth,
        f_partial_ela_avg=f_partial_ela_avg,
        θ_partial_ela_avg=θ_partial_ela_avg,
        d_partial_ela_avg=d_partial_ela_avg,
        de_partial_ela_avg=de_partial_ela_avg,
        d_partial_ela_total_avg=d_partial_ela_total_avg,
        f_ela_avg=f_ela_avg,
        θ_ela_avg=θ_ela_avg,
        d_ela_avg=d_ela_avg,
        de_ela_avg=de_ela_avg,
        d_ela_total_avg=d_ela_total_avg,
        AggrPartialSemiElaInvE_nonM=AggrPartialSemiElaInvE_nonM,
        AggrPartialSemiElaInv_M=AggrPartialSemiElaInv_M,
        AggrPartialSemiElaInvCW=AggrPartialSemiElaInvCW,
        AggrPartialSemiElaInvCE=AggrPartialSemiElaInvCE,
        AggrPartialSemiElaInvE=AggrPartialSemiElaInvE,
        AggrPartialSemiElaInvC_M=AggrPartialSemiElaInvC_M,
        AggrPartialSemiElaInv=AggrPartialSemiElaInv,
        PrivateFirmOwnCapByNWShare=PrivateFirmOwnCapByNWShare,
        PrivateFirmCapByNWShare=PrivateFirmCapByNWShare,
        IlliqAssShare=IlliqAssShare)

    return MITStats

end

# Keep the grouping and economic aggregation in the statistics layer.
function capital_shares_by_net_worth(SSall, PolE, G, P)
    bounds = [-Inf, 2.0, 6.0, 10.0, Inf] ./ P.NomGDPFac .* 1e6
    net_worth_e = reshape(G.b_grid_e, :, 1, 1, 1) .+
        reshape(G.a_grid, 1, :, 1, 1) .+ reshape(G.e_grid, 1, 1, :, 1)
    net_worth_w = reshape(G.b_grid, :, 1, 1, 1) .+
        reshape(G.a_grid, 1, :, 1, 1)
    own_capital_mass = reshape(G.e_grid, 1, 1, :, 1) .* SSall.ge .* P.share_entre
    private_capital_mass = reshape(PolE.θ, size(SSall.ge)) .* own_capital_mass
    worker_capital_mass = reshape(G.a_grid, 1, :, 1, 1) .* SSall.g .* (1 - P.share_entre)
    grouped_share(mass, nw, total) = [sum(mass .* ((bounds[i] .<= nw) .&
        (nw .< bounds[i+1]))) / total for i in 1:length(bounds)-1]
    private_shares = grouped_share(private_capital_mass, net_worth_e,
        SSall.AggrFirmSav + SSall.AggrCapDemE)
    isapprox(sum(private_shares), 1; atol=1e-8) ||
        error("Private capital shares must sum to one")
    return (
        PrivateFirmOwnCapByNWShare=grouped_share(own_capital_mass,
            net_worth_e, SSall.AggrFirmSav),
        PrivateFirmCapByNWShare=private_shares,
        IlliqAssShare=grouped_share(worker_capital_mass,
            net_worth_w, SSall.AggrIlliqSav),
    )
end

function distributionavg(x, g)
    ll = length(g[1, :])
    avg = Vector{Float64}(undef, ll)
    for l=1:ll
        avg[l] = dot(x[:], g[:, l]) / sum(g[:, l])
    end
    return avg
end

function percmask(g, bounds_lo, bounds_up, ValueMat)
    ll = length(bounds_up)
    perc_up = Vector{Float64}(undef, ll)
    perc_lo = Vector{Float64}(undef, ll)
    sortindex = sortperm(vec(ValueMat))
    ValueMatSorted = ValueMat[sortindex]
    g_unsort = g ./ sum(g)
    g_sort = g_unsort[sortindex]
    cdf = cumsum(g_sort)
    cdf = cdf ./ maximum(cdf)
    for gg in eachindex(bounds_up)
        perc_lo[gg] = ValueMatSorted[findfirst(cdf .>= bounds_lo[gg])]
        perc_up[gg] = ValueMatSorted[findfirst(cdf .>= bounds_up[gg])]
    end
    partconsidered = valuemask(g, perc_lo, perc_up, ValueMat)
    return partconsidered, perc_up, perc_lo
end

function valuemask(g, bounds_lo, bounds_up, ValueMat)
    partconsidered = repeat(g[:], 1, length(bounds_up))
    for gg in eachindex(bounds_up)
        partconsidered[(ValueMat[:] .< bounds_lo[gg]) .| (ValueMat[:] .>= bounds_up[gg]), gg] .= 0.0
    end
    return partconsidered
end

function change_bygroups(initial_states, initial_statevar, next_statevar_steady, next_statevar_MIT, next_statevar_Partial, perc_lo, perc_up, ValueMat)

    change_steady = next_statevar_steady .- initial_statevar
    change_MIT = next_statevar_MIT .- initial_statevar
    change_Partial = next_statevar_Partial .- initial_statevar

    effect_MIT = change_MIT - change_steady
    effect_Partial = change_Partial - change_steady

    AvgEffect_MIT = mean(effect_MIT)
    AvgEffect_Partial = mean(effect_Partial)

    initial_value = ValueMat[initial_states]

    NoGroups = length(perc_lo)

    AvgEffect_MIT_bygroup = fill(NaN, NoGroups)
    AvgEffect_Partial_bygroup = fill(NaN, NoGroups)
    Nsim_bygroup = zeros(Int, NoGroups)

    for i = 1:NoGroups
        mask = (initial_value .>= perc_lo[i]) .&
               (initial_value .< perc_up[i])

        Nsim_bygroup[i] = count(mask)

        if Nsim_bygroup[i] > 0
            AvgEffect_MIT_bygroup[i] = mean(effect_MIT[mask])./mean(initial_statevar[mask])
            AvgEffect_Partial_bygroup[i] = mean(effect_Partial[mask])./mean(initial_statevar[mask])
        end
    end

    return AvgEffect_MIT_bygroup, AvgEffect_Partial_bygroup, AvgEffect_MIT, AvgEffect_Partial, Nsim_bygroup

end


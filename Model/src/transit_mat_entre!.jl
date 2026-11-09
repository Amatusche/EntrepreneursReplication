function transit_mat_entre!(Aentre, Ic_B, Ic_F, Id_B, Id_F, If_B, If_F,
    sce_B, sce_F, sde_B, sde_F, sf_B, sf_F, dB, dF, fB, fF,
    adrift_e, edrift, FirmIncomeConstants, P, G, Opt, TempE::EntreTemp)

    Xb = TempE.trans_Xb
    Y = TempE.trans_Y
    Zb = TempE.trans_Zb
    Xa = TempE.trans_Xa
    Za = TempE.trans_Za
    Xe = TempE.trans_Xe
    Ze = TempE.trans_Ze

    # These two state-sized buffers are first used for the a drift and then
    # overwritten for the e drift.
    backward_drift = TempE.trans_B
    forward_drift = TempE.trans_F

    # Liquid-asset transitions
    @views begin
        Xb .= -1.0 .* Ic_B[2:G.Ne] .* sce_B[2:G.Ne]
        Xb .-= Id_B[2:G.Ne] .* sde_B[2:G.Ne]
        Xb .-= If_B[2:G.Ne] .* sf_B[2:G.Ne]
        Xb ./= G.rep_deltaE_b_B[2:G.Ne]

        Y .= (Ic_B[:] .* sce_B[:] .+ Id_B[:] .* sde_B[:] .+ If_B[:] .* sf_B[:]) ./ G.rep_deltaE_b_B .-
             (Ic_F[:] .* sce_F[:] .+ Id_F[:] .* sde_F[:] .+ If_F[:] .* sf_F[:]) ./ G.rep_deltaE_b_F

        Zb .= (Ic_F[1:G.Ne-1] .* sce_F[1:G.Ne-1] .+ Id_F[1:G.Ne-1] .* sde_F[1:G.Ne-1] .+
               If_F[1:G.Ne-1] .* sf_F[1:G.Ne-1]) ./ G.rep_deltaE_b_F[1:G.Ne-1]
    end

    # Passive illiquid-asset transitions
    backward_drift .= min.(dB, 0.0)
    forward_drift .= max.(dF, 0.0) .+ adrift_e
    @views forward_drift[:, end, :, :] .= 0.0

    @views begin
        Xa .= -1.0 .* backward_drift[G.N_be+1:G.Ne] ./ G.rep_deltaE_a_B[G.N_be+1:G.Ne]
        Y .+= backward_drift[:] ./ G.rep_deltaE_a_B .- forward_drift[:] ./ G.rep_deltaE_a_F
        Za .= forward_drift[1:G.Ne-G.N_be] ./ G.rep_deltaE_a_F[1:G.Ne-G.N_be]
    end

    # Business-equity transitions, reusing the same two drift buffers
    expected_depreciation_rate = P.delta_e + P.lambda_e / (1 + P.eta) * Opt.entreriskjump
    nonliquidation_constant = -P.delta_e + P.tauE * expected_depreciation_rate
    # Positive operating income is split according to intoliqshare. Operating
    # losses are borne by business equity rather than discarded by the liquid
    # limited-liability floor. Liquidation income is always paid out in liquid
    # form and therefore does not enter the e drift.
    @. backward_drift = min(fB, 0.0) +
        min(edrift - (FirmIncomeConstants - nonliquidation_constant * G.e_grid_MAT), 0.0)
    @. forward_drift = max(fF, 0.0) + (1.0 - P.intoliqshare) *
        max(edrift - (FirmIncomeConstants - nonliquidation_constant * G.e_grid_MAT), 0.0)
    @views forward_drift[:, :, end, :] .= 0.0

    @views begin
        Xe .= -1.0 .* backward_drift[G.N_abe+1:G.Ne] ./ G.rep_deltaE_e_B[G.N_abe+1:G.Ne]
        Y .+= backward_drift[:] ./ G.rep_deltaE_e_B .- forward_drift[:] ./ G.rep_deltaE_e_F
        Ze .= forward_drift[1:G.Ne-G.N_abe] ./ G.rep_deltaE_e_F[1:G.Ne-G.N_abe]
    end

    Aentreassemble!(Aentre, Y, Xa, Xe, Xb, Za, Ze, Zb, G)

    return nothing
end

"""
    transit_mat_fixed_entrepreneur!(Aentre, PolE, RbE, rk, TransE,
        RealGrossProfitsnoLevPerCap, FirmIncomeConstants, adriftE,
        P, G, Opt, TempE)

Build the entrepreneurial transition matrix from fixed consumption, portfolio,
business-investment, and leverage policies, while evaluating all realized
budget flows at actual prices. This is used for entrepreneurial counterfactuals
whose decisions are formed using perceived steady-state prices.
"""
function transit_mat_fixed_entrepreneur!(Aentre, PolE, RbE, rk, TransE,
    RealGrossProfitsnoLevPerCap, FirmIncomeConstants, adriftE,
    P, G, Opt, TempE::EntreTemp)

    Xb = TempE.trans_Xb
    Y = TempE.trans_Y
    Zb = TempE.trans_Zb
    Xa = TempE.trans_Xa
    Za = TempE.trans_Za
    Xe = TempE.trans_Xe
    Ze = TempE.trans_Ze

    consumption_drift = TempE.sce_B
    deposit_liquid_drift = TempE.sde_B
    investment_liquid_drift = TempE.sf_B
    business_income = TempE.edrift0
    backward_drift = TempE.trans_B
    forward_drift = TempE.trans_F

    ce = reshape(PolE.ce, G.N_be, G.N_a, G.N_e, G.N_y)
    de = reshape(PolE.de, G.N_be, G.N_a, G.N_e, G.N_y)
    f = reshape(PolE.f, G.N_be, G.N_a, G.N_e, G.N_y)
    theta = reshape(PolE.θ, G.N_be, G.N_a, G.N_e, G.N_y)

    leverage_power = (P.nu * P.alphaaE) /
                     (1.0 - P.nu * (1.0 - P.alphaaE))
    @. business_income = FirmIncomeConstants +
        (1.0 - P.tauE) * RealGrossProfitsnoLevPerCap * G.e_grid_MAT *
        theta^leverage_power -
        (1.0 - P.tauE) * (theta - 1.0) * G.e_grid_MAT * rk

    # Liquid-account flows at actual prices. The policies themselves remain
    # those selected in the perceived-price problem.
    scebuild!(consumption_drift, RbE, TransE, ce, business_income,
        FirmIncomeConstants, G, P, Opt)
    sdebuild!(deposit_liquid_drift, de, G.a_grid_MAT_e_adj,
        P.χ0, P.χ1, P.χ2)
    sdebuild!(investment_liquid_drift, f, G.e_grid_MAT_adj,
        0.0, P.χe1, P.χe2)

    @. backward_drift = min(consumption_drift, 0.0) +
        min(deposit_liquid_drift, 0.0) + min(investment_liquid_drift, 0.0)
    @. forward_drift = max(consumption_drift, 0.0) +
        max(deposit_liquid_drift, 0.0) + max(investment_liquid_drift, 0.0)
    @views begin
        backward_drift[1, :, :, :] .= 0.0
        forward_drift[end, :, :, :] .= 0.0

        Xb .= -backward_drift[2:G.Ne] ./ G.rep_deltaE_b_B[2:G.Ne]
        Y .= backward_drift[:] ./ G.rep_deltaE_b_B .-
             forward_drift[:] ./ G.rep_deltaE_b_F
        Zb .= forward_drift[1:G.Ne-1] ./ G.rep_deltaE_b_F[1:G.Ne-1]
    end

    # Illiquid financial wealth.
    @. backward_drift = min(de, 0.0)
    @. forward_drift = max(de, 0.0) + adriftE
    @views begin
        backward_drift[:, 1, :, :] .= 0.0
        forward_drift[:, end, :, :] .= 0.0

        Xa .= -backward_drift[G.N_be+1:G.Ne] ./
              G.rep_deltaE_a_B[G.N_be+1:G.Ne]
        Y .+= backward_drift[:] ./ G.rep_deltaE_a_B .-
              forward_drift[:] ./ G.rep_deltaE_a_F
        Za .= forward_drift[1:G.Ne-G.N_be] ./
              G.rep_deltaE_a_F[1:G.Ne-G.N_be]
    end

    # Private business equity. Liquidation proceeds enter the liquid account;
    # operating losses reduce business equity and retained positive operating
    # income increases it.
    expected_depreciation_rate = P.delta_e +
                                 P.lambda_e / (1.0 + P.eta) * Opt.entreriskjump
    nonliquidation_constant = -P.delta_e +
                              P.tauE * expected_depreciation_rate
    @. backward_drift = min(f, 0.0) +
        min(business_income -
            (FirmIncomeConstants - nonliquidation_constant * G.e_grid_MAT), 0.0)
    @. forward_drift = max(f, 0.0) + (1.0 - P.intoliqshare) *
        max(business_income -
            (FirmIncomeConstants - nonliquidation_constant * G.e_grid_MAT), 0.0)
    @views begin
        backward_drift[:, :, 1, :] .= 0.0
        forward_drift[:, :, end, :] .= 0.0

        Xe .= -backward_drift[G.N_abe+1:G.Ne] ./
              G.rep_deltaE_e_B[G.N_abe+1:G.Ne]
        Y .+= backward_drift[:] ./ G.rep_deltaE_e_B .-
              forward_drift[:] ./ G.rep_deltaE_e_F
        Ze .= forward_drift[1:G.Ne-G.N_abe] ./
              G.rep_deltaE_e_F[1:G.Ne-G.N_abe]
    end

    Aentreassemble!(Aentre, Y, Xa, Xe, Xb, Za, Ze, Zb, G)

    # For risky leverage, diffusion scales with the square of the fixed
    # leverage policy. Risky jumps are added later by complete_entre_generator.
    if Opt.leverage_mode == :risky
        for ii = eachindex(G.iadiff)
            state = G.iadiff[ii]
            pos = G.diff_nzidx[ii]
            Aentre.nzval[pos] += (PolE.θ[state]^2.0 - 1.0) * G.vadiff[ii]
        end
    end

    return nothing
end

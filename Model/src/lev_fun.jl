"""
`lev_fun0(θ, FirmIncomeConstantsxx, RealGrossProfitsEnoLevDampedxx, a, rk, rbe, b, ke, TransExx, P, ve_exx)`

This function evaluates the FOC for the leverage choice θ when the b drift is zero.
"""
function lev_fun0(θ, FirmIncomeConstantsxx, RealGrossProfitsnoLevPerCap, a,
    rk, rbe, b, ke, TransExx, P, ve_eexx, mode,
    bb, aa, ee, yy, ve_e_F_raw, G, Opt)

    expected_depreciation_rate = P.delta_e + P.lambda_e / (1 + P.eta) * Opt.entreriskjump
    nonliquidation_constant = -P.delta_e + P.tauE * expected_depreciation_rate
    liquidation_income = FirmIncomeConstantsxx - nonliquidation_constant * ke
    operating_income = FirmIncomeConstantsxx - liquidation_income +
                       (1.0 - P.tauE) * RealGrossProfitsnoLevPerCap * ke * θ^a -
                       (1.0 - P.tauE) * (θ - 1.0) * ke * rk
    c = max(liquidation_income + P.intoliqshare * max(operating_income, 0.0) +
            rbe * b + TransExx, eps(Float64))
    marginal_profit = (operating_income > 0.0) * c^(-P.sigma) *
                      P.intoliqshare * (1.0 - P.tauE) *
                      (a * RealGrossProfitsnoLevPerCap * θ^(a - 1.0) - rk)

    if mode == :risky
        diffusion_variance = Opt.entreriskdiffusion ? P.sigma_e^2 : 0.0
        jump_foc = risky_jump_foc(θ, bb, aa, ee, yy, ve_e_F_raw, P, G, Opt)
        return marginal_profit + diffusion_variance * θ * ve_eexx + jump_foc
    end
    return marginal_profit
end

"""
`lev_fun(θ, RealGrossProfitsnoLevPerCap, a, rk, P, ve_exx, ve_b)`
This function evaluates the FOC for the leverage choice θ when the b drift is non-zero.
"""
function lev_fun(θ, RealGrossProfitsnoLevPerCap, a, rk, P, ve_eexx, ve_b,
    mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt)

    marginal_profit = ve_b * P.intoliqshare * (1.0 - P.tauE) *
                      (a * RealGrossProfitsnoLevPerCap * θ^(a - 1.0) - rk)

    if mode == :risky
        diffusion_variance = Opt.entreriskdiffusion ? P.sigma_e^2 : 0.0
        jump_foc = risky_jump_foc(θ, bb, aa, ee, yy, ve_e_F_raw, P, G, Opt)
        return marginal_profit + diffusion_variance * θ * ve_eexx + jump_foc
    end
    return marginal_profit
end

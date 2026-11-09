function cebuild!(ce, RbE, Ic_0, b_grid_MAT_e, edrift, FirmIncomeConstants,
    TransE, G, P, Opt)
    ind = CartesianIndices(Ic_0)
    bb = 0
    expected_depreciation_rate = P.delta_e + P.lambda_e / (1 + P.eta) * Opt.entreriskjump
    nonliquidation_constant = -P.delta_e + P.tauE * expected_depreciation_rate
    for xx = findall(Ic_0[:])
    bb = ind[xx][1]
    liquidation_income = FirmIncomeConstants[xx] -
                         nonliquidation_constant * G.e_grid_MAT[xx]
    operating_income = edrift[xx] - liquidation_income
    ce[xx] = RbE[bb] .* b_grid_MAT_e[xx] .+
             liquidation_income + P.intoliqshare .* max(operating_income, 0.0) .+
             TransE[xx]
    if ce[xx] < 0.0
        ce[xx] = Opt.min_c
    end
end
end

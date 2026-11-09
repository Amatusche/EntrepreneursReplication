function scebuild!(sce_B, Rb, TransE, ce_B, edrift, FirmIncomeConstants,
    G, P, Opt)
ind = CartesianIndices(sce_B)
bb = 0
expected_depreciation_rate = P.delta_e + P.lambda_e / (1 + P.eta) * Opt.entreriskjump
nonliquidation_constant = -P.delta_e + P.tauE * expected_depreciation_rate
for xx = 1:G.Ne
    bb = ind[xx][1]
    liquidation_income = FirmIncomeConstants[xx] -
                         nonliquidation_constant * G.e_grid_MAT[xx]
    operating_income = edrift[xx] - liquidation_income
    # Limited liability applies to operating income, but liquidation proceeds
    # are a return of principal and always accrue to the liquid account.
    sce_B[xx] = Rb[bb] .* G.b_grid_MAT_e[xx] .+ TransE[xx] .+
                liquidation_income + P.intoliqshare .* max(operating_income, 0.0) .-
                ce_B[xx]
end
end

"`MITall!(MIT, n, TransPolE, TransPolW, ge, g)`

Computes adjustment costs etc. from MIT structure and policy functions.
"
function MITall!(MIT, n, TransPolE, TransPolW, ge, g, P, G, Opt,
    TransPolWNA=false, gNA=false, shareNA=0.0)

        ge_vec = vec(ge)
        g_vec = vec(g)
        a_e_vec = vec(G.a_grid_MAT_e_adj)
        a_w_vec = vec(G.a_grid_MAT_adj)
        e_vec = vec(G.e_grid_MAT)
        e_adj_vec = vec(G.e_grid_MAT_adj)
        be_vec = vec(G.b_grid_MAT_e)
        b_vec = vec(G.b_grid_MAT)

        MIT.ConsumpE_t[n] = dot(ge_vec, TransPolE.ce) .* P.share_entre
        MIT.ConsumpW_t[n] = dot(g_vec, TransPolW.c) .* (1.0 .- P.share_entre)

        portf_adj_e = 0.0
        cap_adj_e = 0.0
        inv_e_nonmarket = 0.0
        borrowing_e = 0.0
        entre_denominator = 1.0 - (1.0 - P.alphaaE) * P.nu
        labor_power = (1.0 - P.alphaaE) * P.nu
        capital_power = (P.nu - 1.0) / entre_denominator
        leverage_power = P.nu * P.alphaaE / entre_denominator
        expected_depreciation_rate = P.delta_e +
                                     P.lambda_e / (1.0 + P.eta) * Opt.entreriskjump
        nonliquidation_constant = -P.delta_e +
                                  P.tauE * expected_depreciation_rate
        gross_profit_constant = (MIT.pe_t[n] * P.entre_prod *
                                 (P.nu * (1.0 - P.alphaaE) / MIT.w_t[n])^labor_power)^
                                (1.0 / entre_denominator) * entre_denominator
        @inbounds for ii in eachindex(ge_vec)
            de = TransPolE.de[ii]
            f = TransPolE.f[ii]
            mass = ge_vec[ii]
            e = e_vec[ii]
            yy = div(ii - 1, G.N_bae) + 1
            gross_profit_per_capital = e > 0.0 ?
                gross_profit_constant * G.y_grid[yy]^(1.0 / entre_denominator) *
                e^capital_power : 1.0
            leverage = TransPolE.θ[ii]
            operating_income = nonliquidation_constant * e +
                (1.0 - P.tauE) * gross_profit_per_capital * e * leverage^leverage_power -
                (1.0 - P.tauE) * (leverage - 1.0) * e * MIT.rk_t[n]
            portf_adj_e += mass * (P.χ0 * abs(de) + P.χ1 * abs(de / a_e_vec[ii]) ^ P.χ2 * a_e_vec[ii])
            cap_adj_e += mass * P.χe1 * abs(f / e_adj_vec[ii]) ^ P.χe2 * e_adj_vec[ii]
            inv_e_nonmarket += mass * (P.delta_e * e + f + min(operating_income, 0.0))
            if be_vec[ii] < 0.0
                borrowing_e += mass * be_vec[ii]
            end
        end
        MIT.PortfAdjE_t[n] = portf_adj_e .* P.share_entre
        MIT.CapAdjE_t[n] = cap_adj_e .* P.share_entre
        MIT.InvE_nonM_t[n] = inv_e_nonmarket .* P.share_entre
        BorrCostsE = -P.kappa .* P.share_entre .* borrowing_e

        portf_adj_w = 0.0
        borrowing_w = 0.0
        @inbounds for ii in eachindex(g_vec)
            d = TransPolW.d[ii]
            mass = g_vec[ii]
            portf_adj_w += mass * (P.χ0 * abs(d) + P.χ1 * abs(d / a_w_vec[ii]) ^ P.χ2 * a_w_vec[ii])
            if b_vec[ii] < 0.0
                borrowing_w += mass * b_vec[ii]
            end
        end
        BorrCostsW = -P.kappa .* (1.0 .- P.share_entre) .* borrowing_w
        MIT.PortfAdjW_t[n] = portf_adj_w .* (1.0 .- P.share_entre)

        if shareNA > 0.0
            gNA_vec = vec(gNA)
            ConsumpWNA_t = dot(gNA_vec, TransPolWNA.c) .* (1.0 .- P.share_entre)
            MIT.ConsumpW_t[n] = MIT.ConsumpW_t[n]*(1.0.-shareNA) +  ConsumpWNA_t*shareNA
            borrowing_w_na = 0.0
            portf_adj_w_na = 0.0
            @inbounds for ii in eachindex(gNA_vec)
                d = TransPolWNA.d[ii]
                mass = gNA_vec[ii]
                portf_adj_w_na += mass * (P.χ0 * abs(d) + P.χ1 * abs(d / a_w_vec[ii]) ^ P.χ2 * a_w_vec[ii])
                if b_vec[ii] < 0.0
                    borrowing_w_na += mass * b_vec[ii]
                end
            end
            BorrCostsW = BorrCostsW * (1.0 - shareNA) - P.kappa * (1.0 - P.share_entre) * borrowing_w_na * shareNA
            PortfAdjWNA_t = portf_adj_w_na * (1.0 - P.share_entre)
            MIT.PortfAdjW_t[n] = MIT.PortfAdjW_t[n]*(1.0-shareNA) .+ PortfAdjWNA_t*shareNA
        end

        MIT.PortfAdj_t[n] =  MIT.PortfAdjW_t[n].+  MIT.PortfAdjE_t[n]
        MIT.C_t[n] = MIT.ConsumpW_t[n] .+ MIT.ConsumpE_t[n]
        MIT.BorrCosts_t[n] =  BorrCostsW .+ BorrCostsE
        MIT.KW_t[n] = MIT.AggrCapSup_t[n]*MIT.AggrIlliqSavW_t[n]/MIT.AggrIlliqSav_t[n]
        MIT.KE_t[n] = MIT.AggrCapSup_t[n]*MIT.AggrIlliqSavE_t[n]/MIT.AggrIlliqSav_t[n]

    if n > 1
        MIT.InvCE_t[n-1] = (MIT.KE_t[n] .- MIT.KE_t[n-1]) ./ G.dt[n-1] .+ P.delta .* MIT.KE_t[n-1]
        MIT.InvCW_t[n-1] = (MIT.KW_t[n] .- MIT.KW_t[n-1]) ./ G.dt[n-1] .+ P.delta .* MIT.KW_t[n-1]
    end

end

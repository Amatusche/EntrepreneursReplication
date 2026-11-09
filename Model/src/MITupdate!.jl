" `function MITupdate!(MIT, AggrOutputE_new, AggrLabDemE_new, tax_rev_entre_new, g_tilde, g_tilde_e, SSallend, A_t, Aentre_local_t, TransPolW, TransPolE, g, ge, TempW, TempE, P, G, Opt, entre_exp, worker_exp)`

    Performs one updating step of MIT shock. The variable `entre_exp` indicates how entrepreneurs behave.
    `entre_exp = 0`: Baseline model.
    `entre_exp = 1`: all entrepreneurial variables unchanged, as given by `MIT` -> no market clearing because entrepreneurs' budget constraint does not hold
    `entre_exp = 2`: entrepreneurs act as if all prices were in steady state, but incomes still adjust -> market clearing
    `entre_exp = 3`: same as 2 but only liquid rate perceived to be in steady state. -> market clearing

    worker_exp = 4 some workers do not respond to rb change
"
function MITupdate!(MIT, AggrOutputE_new, AggrLabDemE_new, tax_rev_entre_new,
    g_tilde, g_tilde_e, SSallend, A_t, Aentre_local_t, TransPolW, TransPolE, g, ge,
    TempW, TempE,
    P, G, Opt, entre_exp, worker_exp, g_tildeNA=false, ANA_t=false, TransPolWNA=false, gNA=false, shareNA=0.0)

    Bentre = deepcopy(Aentre_local_t[1])
    B = deepcopy(A_t[1])
    minusAT_e = deepcopy(Aentre_local_t[1])
    minusATassemble!(minusAT_e, Aentre_local_t[1], P.zeta, G.iaexo, G.jaexo)
    minusAT = deepcopy(A_t[1])
    minusATassemble!(minusAT, A_t[1], P.zeta, G.iaz, G.jaz)

    v = deepcopy(SSallend.v)     # at beginning of each new round of iterations impose that terminal value function is SS.v from steady state
    if P.share_entre > 0.0
        ve = deepcopy(SSallend.ve)
    end

    Rb = Array{Float64}(undef, G.N_b)
    RbE = Array{Float64}(undef, G.N_be)
    TransW = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    lab0 = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    TransE = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    re = Array{Float64}(undef, G.N_e, G.N_y)
    b = TempW.b
    be = TempE.b
    LabSupBAJK = Array{Float64}(undef, G.N)
    LabDemENoLev = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    LabDemE = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    OutputEnoLev = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    OutputE = deepcopy(SSallend.OutputE)
    RealGrossProfitsnoLevPerCap = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    FirmIncomeConstants = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    RealGrossProfitsE = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    expected_depreciation_rate = P.delta_e + P.lambda_e / (1 + P.eta) * Opt.entreriskjump

    g_marg_b = Array{Float64}(undef, G.N_b)
    g_marg_a = Array{Float64}(undef, G.N_a)
    ge_marg_b = Array{Float64}(undef, G.N_be)
    ge_marg_a = Array{Float64}(undef, G.N_a)
    ge_marg_e = Array{Float64}(undef, G.N_e)

    w = Float64
    what = Float64
    ra = Float64
    rb = Float64
    pe = Float64
    # dampen illiquid returns
    ltau0a = Float64
    adriftE = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    adrift = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    if entre_exp >= 1 || worker_exp >= 1
        RbSSE = SSallend.rb .+ P.kappa .* (G.b_grid_e .< 0.0) .+ P.zeta .* SSallend.zeta_corr
        RbSS = SSallend.rb .+ P.kappa .* (G.b_grid .< 0.0) .+ P.zeta .* SSallend.zeta_corr
        RealGrossProfitsnoLevPerCapSS = (SSallend.pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (P.nu * (1.0 - P.alphaaE) / SSallend.w)^((1 - P.alphaaE) * P.nu)) .^ (1 / (1 - (1 - P.alphaaE) * P.nu)) .* (1 - (1 - P.alphaaE) * P.nu) .* G.e_grid_MAT .^ ((P.nu - 1) / (1 - (1 - P.alphaaE) * P.nu))
        RealGrossProfitsnoLevPerCapSS[G.e_grid_MAT.==0.0] .= 1.0
        ltau0a = SSallend.ra .* (G.a_max .* 0.999) .^ (1.0 .- Opt.ltau)
        adriftSSE = SSallend.ra .* G.a_grid_MAT_e .- ltau0a .* (G.a_grid_MAT_e .^ Opt.ltau)
        adriftSS = SSallend.ra .* G.a_grid_MAT .- ltau0a .* (G.a_grid_MAT .^ Opt.ltau)
    end

    if worker_exp == 4
        lab0NA = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
        vNA = deepcopy(SSallend.v)
    end

    for n = G.N_t:-1:1
        if mod(n / 50, 1) == 0 && Opt.display > 1
            println("   going backward, period " * string(n))
        end
        Rb .= MIT.rb_t[n] .+ P.kappa .* (G.b_grid .< 0.0) .+ P.zeta .* MIT.zeta_corr_t[n]
        RbE .= MIT.rb_t[n] .+ P.kappa .* (G.b_grid_e .< 0.0) .+ P.zeta .* MIT.zeta_corr_t[n]
        FirmIncomeConstants .= (-P.delta_e + P.tauE * expected_depreciation_rate + P.zeta * MIT.zeta_corr_t[n]) .* G.e_grid_MAT

        w = MIT.w_t[n]
        what = MIT.what_t[n]
        ra = MIT.ra_t[n]
        rb = MIT.rb_t[n]
        rk = MIT.rk_t[n]

        MIT.RealProfitsToW_t[n] = (1.0 .- P.omega_liq) .* MIT.RealProfits_t[n] .+ MIT.UnionProfits_t[n] # total profits that go to workers liquid account in proportion to labor income
        MIT.RealProfitsToE_t[n] = P.share_entre > 0.0 ? P.omega_liq .* MIT.omega_entre[n] .* MIT.RealProfits_t[n] : 0.0 # total profits that go to entrepreneurs in proportion to output income
        MIT.AvgPerCapitaRealProfitsToW_t[n] = MIT.RealProfitsToW_t[n] ./ (1.0 .- P.share_entre) # average profits PER WORKER going to liquid account
        MIT.AvgPerCapitaRealProfitsToE_t[n] = P.share_entre > 0.0 ? MIT.RealProfitsToE_t[n] ./ P.share_entre : 0.0 # average profits PER ENTREPRENEURS going to business income

        TransW .= MIT.Trans_t[n] .+ (1.0 .- P.tau) .* MIT.AvgPerCapitaRealProfitsToW_t[n] .* G.z_grid_MAT ./ G.AvgLabProduct
        if P.share_entre > 0.0
            pe = MIT.pe_t[n]
            LabDemENoLev .= (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (1.0 .- P.alphaaE) .* P.nu .* G.e_grid_MAT .^ (P.alphaaE * P.nu) ./ w) .^ (1.0 ./ (1.0 .- (1.0 .- P.alphaaE) .* P.nu))
            OutputEnoLev .= P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (G.e_grid_MAT .^ P.alphaaE .* LabDemENoLev .^ (1.0 .- P.alphaaE)) .^ P.nu
            vec(OutputE) .= vec(OutputEnoLev) .* TransPolE[n].θ .^ ((P.alphaaE .* P.nu) ./ (1.0 .- P.nu .* (1.0 .- P.alphaaE)))
            TransE .= MIT.Trans_t[n] .+ (1.0 .- P.tauE) .* MIT.AvgPerCapitaRealProfitsToE_t[n] .* OutputE ./ (MIT.AggrOutputE_t[n] ./ P.share_entre)
            RealGrossProfitsnoLevPerCap .= (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (P.nu * (1.0 - P.alphaaE) / w)^((1 - P.alphaaE) * P.nu)) .^ (1 / (1 - (1 - P.alphaaE) * P.nu)) .* (1 - (1 - P.alphaaE) * P.nu) .* G.e_grid_MAT .^ ((P.nu - 1) / (1 - (1 - P.alphaaE) * P.nu))
            RealGrossProfitsnoLevPerCap[G.e_grid_MAT .== 0.0] .= 1.0
        end

        # dampen illiquid returns
        ltau0a = ra .* (G.a_max .* 0.999) .^ (1.0 .- Opt.ltau)
        adriftE .= ra .* G.a_grid_MAT_e .- ltau0a .* (G.a_grid_MAT_e .^ Opt.ltau)
        adrift .= ra .* G.a_grid_MAT .- ltau0a .* (G.a_grid_MAT .^ Opt.ltau)

        if worker_exp == 0
            lab0update!(lab0, TransW, what, Rb, G, P, Opt)
            HHsolve!(v, TransPolW[n], A_t[n], B, Rb, what, TransW, lab0, adrift,
                P, G, Opt, G.dt[n], TempW)
        elseif worker_exp == 3 #  workers do not respond to Rb change
            # In this case, we use the steady state rb to compute decisions -> update value function because other prices change
            lab0update!(lab0, TransW, what, RbSS, G, P, Opt)
            HHsolve!(v, TransPolW[n], A_t[n], B, RbSS, what, TransW, lab0, adrift,
                P, G, Opt, G.dt[n], TempW)
            # Keep the policies chosen under the steady-state liquid return,
            # but evolve the distribution using budget drifts at the actual

            # zero-drift indicators in the realized transition matrix.
            TransPolW[n].d .= vec(TempW.Id_B) .* vec(TempW.d_B) .+
                              vec(TempW.Id_F) .* vec(TempW.d_F)
            transit_mat_fixed_worker!(A_t[n], TransPolW[n], Rb, what, TransW, adrift, P, G, TempW)
        elseif worker_exp == 4 # some workers do no respond to rb change
            # In this case, we use the steady state rb to compute decisions -> update value function because other prices change
            lab0update!(lab0NA, TransW, what, RbSS, G, P, Opt)
            HHsolve!(vNA, TransPolWNA[n], ANA_t[n], B, RbSS, what, TransW, lab0NA, adrift,
                P, G, Opt, G.dt[n], TempW)
            # Keep this subgroup's perceived-return policies fixed, but use
            # actual-return budget drifts to evolve its distribution.
            TransPolWNA[n].d .= vec(TempW.Id_B) .* vec(TempW.d_B) .+
                                vec(TempW.Id_F) .* vec(TempW.d_F)
            transit_mat_fixed_worker!(ANA_t[n], TransPolWNA[n], Rb, what, TransW,
                adrift, P, G, TempW)
            # Other workers do respond
            lab0update!(lab0, TransW, what, Rb, G, P, Opt)
            HHsolve!(v, TransPolW[n], A_t[n], B, Rb, what, TransW, lab0, adrift,
                P, G, Opt, G.dt[n], TempW)
        end

        TransPolW[n].d .= vec(TempW.Id_B) .* vec(TempW.d_B) .+ vec(TempW.Id_F) .* vec(TempW.d_F)

        if P.share_entre > 0
            if entre_exp == 0
                HHsolve_entre!(ve, TransPolE[n], Aentre_local_t[n], Bentre, RbE,
                    rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants, adriftE, P, G, Opt, G.dt[n], true, true, TempE)
                TransPolE[n].f .= vec(TempE.If_B) .* vec(TempE.f_B) .+ vec(TempE.If_F) .* vec(TempE.f_F)
                TransPolE[n].de .= vec(TempE.Id_B) .* vec(TempE.d_B) .+ vec(TempE.Id_F) .* vec(TempE.d_F)
            elseif (entre_exp == 1) & (MIT.it[1] == 0)
                # In this case, we use the steady state prices -> no need to update value function
                updatevalfunc = false
                updatedecisions = true
                HHsolve_entre!(ve, TransPolE[n], Aentre_local_t[n], Bentre, RbSSE,
                    SSallend.rk, SSallend.TransE, RealGrossProfitsnoLevPerCapSS, FirmIncomeConstants, adriftSSE, P, G, Opt, G.dt[n], updatevalfunc, updatedecisions, TempE)
                TransPolE[n].f .= vec(TempE.If_B) .* vec(TempE.f_B) .+ vec(TempE.If_F) .* vec(TempE.f_F)
                TransPolE[n].de .= vec(TempE.Id_B) .* vec(TempE.d_B) .+ vec(TempE.Id_F) .* vec(TempE.d_F)
            elseif entre_exp == 2
                # In this case, we use the steady state prices to compute decisions -> no need to update value function
                if MIT.it[1] == 0
                    updatevalfunc = false
                    updatedecisions = true

                    HHsolve_entre!(ve, TransPolE[n], Aentre_local_t[n], Bentre, RbSSE,
                        SSallend.rk, SSallend.TransE, RealGrossProfitsnoLevPerCapSS, FirmIncomeConstants, adriftSSE, P, G, Opt, G.dt[n], updatevalfunc, updatedecisions, TempE)
                    # We base the transitions on actual prices but keep the decisions obtained under the steady state prices
                    TransPolE[n].f .= vec(TempE.If_B) .* vec(TempE.f_B) .+ vec(TempE.If_F) .* vec(TempE.f_F)
                    TransPolE[n].de .= vec(TempE.Id_B) .* vec(TempE.d_B) .+ vec(TempE.Id_F) .* vec(TempE.d_F)
                end
                transit_mat_fixed_entrepreneur!(Aentre_local_t[n], TransPolE[n],
                    RbE, rk, TransE, RealGrossProfitsnoLevPerCap,
                    FirmIncomeConstants, adriftE, P, G, Opt, TempE)
            elseif entre_exp == 3
                # In this case, we use the steady state rb to compute decisions -> update value function because other prices change
                updatevalfunc = true
                updatedecisions = true
                HHsolve_entre!(ve, TransPolE[n], Aentre_local_t[n], Bentre, RbSSE,
                    rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants, adriftE, P, G, Opt, G.dt[n], updatevalfunc, updatedecisions, TempE)
                TransPolE[n].f .= vec(TempE.If_B) .* vec(TempE.f_B) .+ vec(TempE.If_F) .* vec(TempE.f_F)
                TransPolE[n].de .= vec(TempE.Id_B) .* vec(TempE.d_B) .+ vec(TempE.Id_F) .* vec(TempE.d_F)
                transit_mat_fixed_entrepreneur!(Aentre_local_t[n], TransPolE[n],
                    RbE, rk, TransE, RealGrossProfitsnoLevPerCap,
                    FirmIncomeConstants, adriftE, P, G, Opt, TempE)
            end
        end
    end

    for n = 1:G.N_t
        if (mod(n / 50, 1) == 0 || n == 1) && Opt.display > 1
            println("   going forward, period " * string(n))
        end

        # Worker
        b .= g_tilde ./ G.dt[n] .+ G.b_vec .* P.zeta

        x = P.zeta .+ 1.0 ./ G.dt[n]
        minusATassemble!(minusAT, A_t[n], x, G.iaz, G.jaz)
        linear_solve!(g_tilde, minusAT, b, TempW,
            Opt.worker_KFE_solver, Opt.bicgstab_tol;
            kfe=true, label="worker transition KFE")

        g_sum = sum(g_tilde)
        g_tilde .= g_tilde ./ g_sum
        g .= reshape(g_tilde, G.N_b, G.N_a, G.N_j, G.N_k)

        # marginal densities
        sum!(reshape(g_marg_b, G.N_b, 1, 1, 1), g)
        sum!(reshape(g_marg_a, 1, G.N_a, 1, 1), g)

        # compute aggregate savings and labor supply
        MIT.AggrLiqSavW_t[n+1] = dot(g_marg_b, G.b_grid) .* (1.0 .- P.share_entre)
        MIT.AggrIlliqSavW_t[n+1] = dot(g_marg_a, G.a_grid) .* (1.0 .- P.share_entre)

        for jj = 1:G.N
            LabSupBAJK[jj] = G.z_grid_MAT[jj] .* TransPolW[n].lab[jj] .* g[jj]
        end
        MIT.AggrLabSup_t[n] = sum(LabSupBAJK) .* (1.0 .- P.share_entre)

        if worker_exp == 4 # also compute distribution for workers that do not respond, then aggregate
            b .= g_tildeNA ./ G.dt[n] .+ G.b_vec .* P.zeta
            minusATassemble!(minusAT, ANA_t[n], x, G.iaz, G.jaz)
            linear_solve!(g_tildeNA, minusAT, b, TempW,
                Opt.worker_KFE_solver, Opt.bicgstab_tol;
                kfe=true, label="nonresponding-worker transition KFE")

            g_sumNA = sum(g_tildeNA)
            g_tildeNA .= g_tildeNA ./ g_sumNA

            gNA .= reshape(g_tildeNA, G.N_b, G.N_a, G.N_j, G.N_k)

            # marginal densities
            sum!(reshape(g_marg_b, G.N_b, 1, 1, 1), gNA)
            sum!(reshape(g_marg_a, 1, G.N_a, 1, 1), gNA)

            # compute aggregate savings and labor supply
            MIT.AggrLiqSavW_t[n+1] = MIT.AggrLiqSavW_t[n+1] * (1.0 - shareNA) + dot(g_marg_b, G.b_grid) .* (1.0 .- P.share_entre) * shareNA
            MIT.AggrIlliqSavW_t[n+1] = MIT.AggrIlliqSavW_t[n+1] * (1.0 - shareNA) + dot(g_marg_a, G.a_grid) .* (1.0 .- P.share_entre) * shareNA

            for jj = 1:G.N
                LabSupBAJK[jj] = G.z_grid_MAT[jj] .* TransPolWNA[n].lab[jj] .* gNA[jj]
            end
            MIT.AggrLabSup_t[n] = MIT.AggrLabSup_t[n] * (1.0 - shareNA) + sum(LabSupBAJK) .* (1.0 .- P.share_entre) * shareNA

        end

        # Entrepreneurs
        if P.share_entre > 0
            be .= g_tilde_e ./ G.dt[n] .+ G.b_vecE .* P.zeta

            if risky_jump_active(Opt, P, G)
                Aentre_KFE = complete_entre_generator(Aentre_local_t[n], TransPolE[n].θ, P, G, Opt)
                minusAT_e = spdiagm(0 => fill(x, G.Ne)) - sparse(transpose(Aentre_KFE))
            else
                minusATassemble!(minusAT_e, Aentre_local_t[n], x, G.iaexo, G.jaexo)
            end
            rebuild_KFE_preconditioner = risky_jump_active(Opt, P, G)
            linear_solve!(g_tilde_e, minusAT_e, be, TempE,
                Opt.entrepreneur_KFE_solver, Opt.bicgstab_tol;
                kfe=true, rebuild_preconditioner=rebuild_KFE_preconditioner,
                label="entrepreneur transition KFE")

            g_sum = sum(g_tilde_e)
            g_tilde_e .= g_tilde_e ./ g_sum
            ge .= reshape(g_tilde_e, G.N_be, G.N_a, G.N_e, G.N_y)

            # marginal densities
            sum!(reshape(ge_marg_b, G.N_be, 1, 1, 1), ge)
            sum!(reshape(ge_marg_a, 1, G.N_a, 1, 1), ge)
            sum!(reshape(ge_marg_e, 1, 1, G.N_e, 1), ge)

            pe = MIT.pe_t[n]
            w = MIT.w_t[n]
            TransE .= MIT.Trans_t[n] .+ (1.0 .- P.tauE) .* MIT.AvgPerCapitaRealProfitsToE_t[n] .* OutputE ./ (MIT.AggrOutputE_t[n] ./ P.share_entre)
            LabDemENoLev .= (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (1.0 .- P.alphaaE) .* P.nu .* G.e_grid_MAT .^ (P.alphaaE * P.nu) ./ w) .^ (1.0 ./ (1.0 .- (1.0 .- P.alphaaE) .* P.nu))
            OutputEnoLev .= P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (G.e_grid_MAT .^ P.alphaaE .* LabDemENoLev .^ (1.0 .- P.alphaaE)) .^ P.nu

            vec(LabDemE) .= vec(LabDemENoLev) .* TransPolE[n].θ .^ ((P.alphaaE .* P.nu) / (1.0 .- P.nu * (1.0 .- P.alphaaE)))
            vec(OutputE) .= vec(OutputEnoLev) .* TransPolE[n].θ .^ ((P.alphaaE .* P.nu) ./ (1.0 .- P.nu .* (1.0 .- P.alphaaE)))

            AggrOutputE_new[n] = dot(OutputE, ge) .* P.share_entre
            AggrLabDemE_new[n] = dot(LabDemE, ge) .* P.share_entre
            MIT.AggrFirmSav_t[n+1] = dot(G.e_grid, ge_marg_e) .* P.share_entre

            cap_demand_e = 0.0
            @inbounds @simd for ii in eachindex(TransPolE[n].θ, G.e_grid_MAT, ge)
                cap_demand_e += (TransPolE[n].θ[ii] - 1.0) * G.e_grid_MAT[ii] * ge[ii]
            end
            MIT.AggrCapDemE_t[n] = cap_demand_e .* P.share_entre
            # compute aggregate savings
            MIT.AggrLiqSavE_t[n+1] = dot(ge_marg_b, G.b_grid_e) .* P.share_entre
            MIT.AggrIlliqSavE_t[n+1] = dot(ge_marg_a, G.a_grid) .* P.share_entre

            if P.share_entre > 0.0
                vec(RealGrossProfitsE) .= MIT.pe_t[n] .* vec(OutputE) .- MIT.w_t[n] .* vec(LabDemE) .- MIT.rk_t[n] .* vec(G.e_grid_MAT) .* (TransPolE[n].θ .- 1.0)
                taxable_profit = dot(RealGrossProfitsE, ge) - expected_depreciation_rate * dot(G.e_grid_MAT, ge)
                tax_rev_entre_new[n] = P.share_entre .* P.tauE .* taxable_profit .+ P.tauE .* MIT.RealProfitsToE_t[n]
            else
                tax_rev_entre_new[n] = 0.0
            end

        end

        if worker_exp == 4
            MITall!(MIT, n, TransPolE[n], TransPolW[n], ge, g, P, G, Opt,
                TransPolWNA[n], gNA, shareNA)
        else
            MITall!(MIT, n, TransPolE[n], TransPolW[n], ge, g, P, G, Opt)
        end
    end
end

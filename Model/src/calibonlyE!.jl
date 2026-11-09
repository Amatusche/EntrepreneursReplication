function calibonlyE!(SSin, P, G, Opt)
    @unpack v, ve, lab0, TransW, TransE, gtemp, gtempE, rb, ra, pb, pe, pp, w, what, Output_total, AggrOutputE, RealProfits,
    AggrCapSup, CapLabRatio, omega_entre, zeta_corr = SSin

    TempE = EntreTemp(G, Opt)
    ge_tilde = zeros(G.Ne)

    Aentre_local, Bentre = sparse_value_path(G.Aexo, 2)
    Aentre = deepcopy(Aentre_local)
    minusAT_e = deepcopy(Aentre_local)
    minusATassemble!(minusAT_e, Aentre_local, P.zeta, G.iaexo, G.jaexo)

    ce = zeros(G.Ne)
    f = zeros(G.Ne)
    de = zeros(G.Ne)
    θ = ones(G.Ne)
    PolE = PolEstruct(ce, f, de, θ)
    gtempE = zeros(G.Ne)
    RbE = rb .+ P.kappa * (G.b_grid_e .< 0.0) .+ P.zeta * zeta_corr
    mc = 1.0 .- 1.0 ./ P.epsilon + P.theta * P.ρ * P.infl / P.epsilon
    rk = P.alphaa .* mc .* CapLabRatio .^ (P.alphaa .- 1.0)
    expected_depreciation_rate = P.delta_e + P.lambda_e / (1 + P.eta) * Opt.entreriskjump
    FirmIncomeConstants = (-P.delta_e + P.tauE * expected_depreciation_rate + P.zeta * zeta_corr) .* G.e_grid_MAT

    RealGrossProfitsnoLevPerCap = (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (P.nu * (1.0 - P.alphaaE) / w)^((1 - P.alphaaE) * P.nu)) .^ (1 / (1 - (1 - P.alphaaE) * P.nu)) * (1 - (1 - P.alphaaE) * P.nu) .* G.e_grid_MAT .^ ((P.nu - 1) / (1 - (1 - P.alphaaE) * P.nu))
    RealGrossProfitsnoLevPerCap[G.e_grid_MAT.==0.0] .= 1.0

    ltau0a = ra .* (G.a_max .* 0.999) .^ (1 - Opt.ltau)
    adriftE = ra .* G.a_grid_MAT_e .- ltau0a .* (G.a_grid_MAT_e .^ Opt.ltau)
    FirmIncomeConstants = (-P.delta_e + P.tauE * expected_depreciation_rate + P.zeta * zeta_corr) .* G.e_grid_MAT

    iterHJBentre!(ve, PolE, Aentre_local, Bentre, RbE, rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants,
        adriftE, P, G, Opt, TempE)
    PolE.f .= TempE.If_B[:] .* TempE.f_B[:] .+ TempE.If_F[:] .* TempE.f_F[:]
    PolE.de .= TempE.Id_B[:] .* TempE.d_B[:] .+ TempE.Id_F[:] .* TempE.d_F[:]
    Aentre = complete_entre_generator(Aentre_local, PolE.θ, P, G, Opt)
    if risky_jump_active(Opt, P, G)
        minusAT_e = spdiagm(0 => fill(P.zeta, G.Ne)) - sparse(transpose(Aentre))
    else
        minusATassemble!(minusAT_e, Aentre, P.zeta, G.iaexo, G.jaexo)
    end
    rebuild_KFE_preconditioner = risky_jump_active(Opt, P, G)
    linear_solve!(gtempE, minusAT_e, G.b_vecE, TempE,
        Opt.entrepreneur_KFE_solver, Opt.bicgstab_tol;
        kfe=true, rebuild_preconditioner=rebuild_KFE_preconditioner,
        label="entrepreneur KFE")
    ge_tilde .= gtempE ./ sum(gtempE)
    ge = reshape(ge_tilde, G.N_be, G.N_a, G.N_e, G.N_y)

    ge_marg_a = reshape(sum(ge, dims=[4, 3, 1]), G.N_a)
    ge_marg_b = sum(ge, dims=[4, 3, 2])[:]
    ge_marg_e = reshape(sum(ge, dims=[4, 2, 1]), G.N_e)
    ge_marg_bae = reshape(sum(ge, dims=[4]), G.N_be, G.N_a, G.N_e)
    ge_marg_ba = reshape(sum(ge, dims=[4, 3]), G.N_be, G.N_a)
    ge_marg_ey = reshape(sum(ge, dims=[2, 1]), G.N_e, G.N_y)
    AggrIlliqSavE = dot(ge_marg_a, G.a_grid) .* P.share_entre

    Trans = P.transshare .* Output_total
    borrowing_costs_worker_MAT = -P.kappa .* G.b_grid_MAT .* (G.b_grid_MAT .< 0.0)
    borrowing_costs_entre_MAT = -P.kappa .* G.b_grid_MAT_e .* (G.b_grid_MAT_e .< 0.0)
    adj_costs_a_worker_MAT = P.χ0 .* abs.(Pol.d) .+ P.χ1 .* abs.(Pol.d ./ G.a_grid_MAT_adj[:]) .^ P.χ2 .* G.a_grid_MAT_adj[:]

    adj_costs_entre_e_E = P.χe1 .* abs.(PolE.f ./ G.e_grid_MAT_adj[:]) .^ P.χe2 .* G.e_grid_MAT_adj[:]
    adj_costs_entre_a_E = P.χ0 .* abs.(PolE.de) .+ P.χ1 .* abs.(PolE.de ./ G.a_grid_MAT_e_adj[:]) .^ P.χ2 .* G.a_grid_MAT_e_adj[:]
    PortfAdjE = P.share_entre .* dot(ge[:], adj_costs_entre_a_E[:])
    CapAdjE = P.share_entre .* dot(ge[:], adj_costs_entre_e_E[:])
    BorrCostsE = P.share_entre .* dot(ge[:], borrowing_costs_entre_MAT[:])

    LabDemENoLev = (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (1.0 .- P.alphaaE) .* P.nu .* G.e_grid_MAT .^ (P.alphaaE * P.nu) / w) .^ (1.0 ./ (1.0 .- (1.0 .- P.alphaaE) .* P.nu)) # optimal labor demand of entrepreneurs
    LabDemE = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    LabDemE[:] .= LabDemENoLev[:] .* PolE.θ .^ ((P.alphaaE .* P.nu) / (1.0 .- P.nu * (1.0 .- P.alphaaE))) # labor demand of entrepreneurs

    AggrCapDemE = dot((PolE.θ .- 1.0) .* G.e_grid_MAT[:], ge[:]) .* P.share_entre

    OutputEnoLev = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    OutputEnoLev .= P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (G.e_grid_MAT .^ P.alphaaE .* LabDemENoLev .^ (1.0 .- P.alphaaE)) .^ P.nu
    OutputE = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    OutputE[:] .= OutputEnoLev[:] .* PolE.θ .^ ((P.alphaaE .* P.nu) ./ (1.0 .- P.nu .* (1.0 .- P.alphaaE))) # output of entrepreneurs
    AggrOutputE = dot(OutputE[:], ge[:]) .* P.share_entre

    RealGrossProfitsE = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    RealGrossProfitsE[:] .= pe .* OutputE[:] .- w .* LabDemE[:] .- rk .* G.e_grid_MAT[:] .* (PolE.θ .- 1.0)
    AdjCosts = P.theta ./ 2.0 .* P.infl .^ 2.0 .* Output_total
    RealProfits = (1.0 .- pb) .* Output_total .- AdjCosts

    AggrCapDemC = AggrCapSup .- AggrCapDemE
    AggrLabDemC = CapLabRatio .* AggrCapDemC
    AggrLabDemE = dot(LabDemE[:], ge[:]) .* P.share_entre
    AggrLabSup =   AggrLabDemC + AggrLabDemE
    subs_payment_union = P.tauW * w * AggrLabSup # subsidy paid to union
    Tunion = subs_payment_union  # lump-sum tax on union
    UnionProfits = ((1.0 + P.tauW) * w - what) * AggrLabSup - P.thetaw / 2.0 * P.inflw^2.0 * AggrLabSup - Tunion # Remaining profits of union
    RealProfitsToW = (1.0 .- P.omega_liq) .* RealProfits + UnionProfits # total profits that go to workers liquid account in proportion to labor income
    RealProfitsToE = P.omega_liq .* omega_entre .* RealProfits # total profits that go to entrepreneurs in proportion to output income
    AvgPerCapitaRealProfitsToW = RealProfitsToW ./ (1.0 .- P.share_entre) # average profits PER WORKER going to liquid account
    AvgPerCapitaRealProfitsToE = RealProfitsToE ./ P.share_entre # average profits PER ENTREPRENEURS going to business income
    TaxableBusinessIncomeE = RealGrossProfitsE .+ AvgPerCapitaRealProfitsToE .* OutputE ./ (AggrOutputE ./ P.share_entre) .- expected_depreciation_rate .* G.e_grid_MAT

    AggrLiqSavE = dot(ge_marg_b, G.b_grid_e) .* P.share_entre
    AggrFirmSav = dot(ge_marg_e, G.e_grid) .* P.share_entre

    OutputC = AggrLabDemC .^ (1.0 .- P.alphaa) .* AggrCapDemC .^ P.alphaa

    q = RealProfits * (1.0 - omega_entre) * P.omega_liq ./ (rk .- P.delta)

    tax_rev_entre = P.share_entre .* P.tauE .* dot(TaxableBusinessIncomeE[:], ge[:])
    InvC_M = P.delta .* AggrCapDemC
    InvE_M = P.delta .* AggrCapDemE
    Inv_M = InvC_M .+ InvE_M
    OperatingIncomeE = TempE.edrift0 .-
                       P.zeta .* zeta_corr .* G.e_grid_MAT
    OperatingLossE = min.(OperatingIncomeE, 0.0)
    InvE_nonM = P.delta_e .* AggrFirmSav .+
                P.share_entre .* dot(PolE.f .+ OperatingLossE[:], ge[:])
    InvE = InvE_M .+ InvE_nonM
    Inv = Inv_M .+ InvE_nonM

    ConsumpE = dot(ge[:], PolE.ce[:]) .* P.share_entre

    g = zeros(G.N_b, G.N_a, G.N_j, G.N_k)
    g_marg_b = zeros(G.N_b)
    g_marg_a = zeros(G.N_a)
    g_marg_ba = zeros(G.N_b, G.N_a)
    AggrLiqSavW = NaN
    AggrIlliqSavW = NaN
    KW = NaN
    KE = NaN
    InvCE = P.delta * KE
    InvCW = P.delta * KW
    PortfAdjW = NaN
    BorrCostsW = NaN
    BorrCosts = BorrCostsW .+ BorrCostsE
    PortfAdj = PortfAdjW .+ PortfAdjE
    AggrLiqSav = AggrLiqSavW .+ AggrLiqSavE
    AggrIlliqSav = AggrIlliqSavW .+ AggrIlliqSavE
    tax_rev_worker = NaN
    BondDem = AggrLiqSav
    Gexp = tax_rev_worker .+ tax_rev_entre .- rb .* BondDem .- Trans + Tunion - subs_payment_union # government expenditure
    ConsumpW = NaN
    C = ConsumpW .+ ConsumpE

    SSout = SSstruct(v, ve, lab0, TransW, TransE, Trans, gtemp, gtempE, rb, ra, pb, pe, pp, w, what, Output_total, AggrOutputE, RealProfits,
        AggrCapSup, KW, KE, CapLabRatio, omega_entre, zeta_corr, LabDemE, AggrCapDemE, AggrCapDemC, OutputE,
        RealGrossProfitsE, TaxableBusinessIncomeE, g, ge, ge_marg_b, ge_marg_a, ge_marg_e, ge_marg_ey, ge_marg_ba,
        ge_marg_bae, g_marg_b, g_marg_a, g_marg_ba, AggrLiqSavW, AggrIlliqSavW,
        AggrLiqSavE, AggrIlliqSavE, AggrFirmSav, AggrLabDemE, AggrLabSup, AggrLiqSav,
        AggrIlliqSav, AggrLabDemC, rk, OutputC, AdjCosts, PortfAdjW,
        PortfAdjE, CapAdjE, BorrCostsW, BorrCostsE, BorrCosts, PortfAdj,
        BondDem, q, tax_rev_entre, tax_rev_worker, subs_payment_union, Gexp, Inv, Inv_M, InvC_M, InvE_M, InvE_nonM, InvCW, InvCE, InvE, ConsumpE,
        ConsumpW, C, RealProfitsToW, RealProfitsToE, AvgPerCapitaRealProfitsToW, AvgPerCapitaRealProfitsToE, UnionProfits, Tunion)

    SSout, Aentre, PolE
end

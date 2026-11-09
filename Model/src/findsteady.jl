function findsteady(SSInitGuess, P, G, Opt)

    println("Computing steady state")

    if Opt.leverage_mode != :none && (P.intoliqshare < 1.0)
        error("Leverage only implemented for P.intoliqshare = 1.0 (firm profits go into liquid account)")
    end

    @unpack v, ve, lab0, TransW, TransE, gtemp, gtempE, rb, ra, pb, pe, pp, w, what, Output_total, AggrOutputE, RealProfits,
    AggrCapSup, CapLabRatio, omega_entre, zeta_corr = SSInitGuess

    ve = copy(SSInitGuess.ve)

    # Initialize
    TempE = EntreTemp(G, Opt)
    TempW = Temp(G, Opt)

    Aentre_local, Bentre = sparse_value_path(G.Aexo, 2) # local HJB generator and HJB system share their CSC pattern
    Aentre = deepcopy(Aentre_local) # full HJB generator for entrepreneurs including risky jumps
    minusAT_e = deepcopy(Aentre_local)
    minusATassemble!(minusAT_e, Aentre_local, P.zeta, G.iaexo, G.jaexo)

    c = zeros(G.N)
    lab = zeros(G.N)
    ce = zeros(G.Ne)
    f = zeros(G.Ne)
    d = zeros(G.N)
    de = zeros(G.Ne)
    Pol = Polstruct(c, lab, d)
    θ = ones(G.Ne)
    PolE = PolEstruct(ce, f, de, θ)

    A, B = sparse_value_path(G.Az, 2)
    minusAT = deepcopy(A)

    ge_tilde = zeros(G.Ne)
    ge = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    ge_marg_a = Array{Float64}(undef, G.N_a)
    ge_marg_ey = Array{Float64}(undef, G.N_e, G.N_y)

    g_tilde = zeros(G.N)
    g = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    g_marg_a = Array{Float64}(undef, G.N_a)
    ge_marg_b = Array{Float64}(undef, G.N_be)
    g_marg_b = Array{Float64}(undef, G.N_b)
    ge_marg_e = Array{Float64}(undef, G.N_e)
    ge_marg_ba = Array{Float64}(undef, G.N_be, G.N_a)
    g_marg_ba = Array{Float64}(undef, G.N_b, G.N_a)
    ge_marg_bae = Array{Float64}(undef, G.N_be, G.N_a, G.N_e)

    LabSupBAJK = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)

    LabDemENoLev = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    LabDemE = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    OutputEnoLev = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    OutputE = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    RealGrossProfitsnoLevPerCap = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    FirmIncomeConstants = zeros(G.N_be, G.N_a, G.N_e, G.N_y)
    RealGrossProfitsE = zeros(G.N_be, G.N_a, G.N_e, G.N_y)

    Trans = P.transshare * Output_total
    mc = 1.0 .- 1.0 ./ P.epsilon + P.theta * P.ρ * P.infl / P.epsilon
    rk = P.alphaa .* mc .* CapLabRatio .^ (P.alphaa .- 1.0)
    CapLabRatioNew = 0.0

    Rb = rb .+ P.kappa * (G.b_grid .< 0.0) .+ P.zeta * zeta_corr
    RbE = rb .+ P.kappa * (G.b_grid_e .< 0.0) .+ P.zeta * zeta_corr
    lab0 = 1.2 .* ones(G.N_b, G.N_a, G.N_j, G.N_k)

    RealProfitsToW = (1.0 .- P.omega_liq) .* RealProfits # total profits that go to workers liquid account in proportion to labor income
    RealProfitsToE = P.omega_liq .* omega_entre .* RealProfits # total profits that go to entrepreneurs in proportion to output income
    AvgPerCapitaRealProfitsToW = RealProfitsToW ./ (1.0 .- P.share_entre) # average profits PER WORKER going to liquid account
    AvgPerCapitaRealProfitsToE = RealProfitsToE ./ P.share_entre # average profits PER ENTREPRENEURS going to business income

    AggrLiqSav = 0.0
    AggrLabSup = 0.0
    AggrFirmSav = 0.0
    AggrIlliqSavW = 0.0
    AggrIlliqSavE = 0.0
    AggrIlliqSav = 0.0
    AggrLabDemE = 0.0
    AggrLabDemC = 0.0
    AggrCapDemC = 0.0
    OutputC = 0.0
    AdjCosts = 0.0

    expected_depreciation_rate = P.delta_e + P.lambda_e/(1+P.eta).*Opt.entreriskjump

    # Compute optimal capital and labor of entrepreneurs for zero leverage
    if P.share_entre > 0.0
        LabDemENoLev .= (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (1.0 .- P.alphaaE) .* P.nu .* G.e_grid_MAT .^ (P.alphaaE * P.nu) / w) .^ (1.0 ./ (1.0 .- (1.0 .- P.alphaaE) .* P.nu)) # optimal labor demand of entrepreneurs
        OutputEnoLev .= P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (G.e_grid_MAT .^ P.alphaaE .* LabDemENoLev .^ (1.0 .- P.alphaaE)) .^ P.nu

        RealGrossProfitsnoLevPerCap = (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (P.nu * (1.0 - P.alphaaE) / w)^((1 - P.alphaaE) * P.nu)) .^ (1 / (1 - (1 - P.alphaaE) * P.nu)) * (1 - (1 - P.alphaaE) * P.nu) .* G.e_grid_MAT .^ ((P.nu - 1) / (1 - (1 - P.alphaaE) * P.nu))
        RealGrossProfitsnoLevPerCap[G.e_grid_MAT.==0.0] .= 1.0
        AggrCapDemE = 0.0
    else
        LabDemENoLev = zeros(G.N_e, G.N_y)
        AggrCapDemE = 0.0
        RealGrossProfitsnoLevPerCap = zeros(G.N_e, G.N_y)
        OutputEnoLev = zeros(G.N_e, G.N_y)
        ge_marg_ey = zeros(G.N_e, G.N_y)
    end

    # dampen illiquid returns
    ltau0a = ra .* (G.a_max .* 0.999) .^ (1 - Opt.ltau)
    adriftE = ra .* G.a_grid_MAT_e .- ltau0a .* (G.a_grid_MAT_e .^ Opt.ltau)
    adrift = ra .* G.a_grid_MAT .- ltau0a .* (G.a_grid_MAT .^ Opt.ltau)
    if P.share_entre > 0.0
        FirmIncomeConstants .=(-P.delta_e + P.tauE * expected_depreciation_rate + P.zeta * zeta_corr) .* G.e_grid_MAT # Constant firm income terms: tax deduction for depreciation and annuity (from accidental bequests)
    end
    # P.tauE * expected_depreciation_rate captures the tax deductions for average depreciation. Only the deterministic part of the depreciation is reinvested automatically and reduces the cash flow from the firm.

    # Description of algorithm:
    # Start with a guess for interest rate and labor (demand).
    # 1. Compute wage from interest rate and marginal costs (firm problem, Phillips curve)
    # 2. Compute capital demand and output and use these to compute transfers
    #       etc. that are relevant for household problem.
    # 2. Solve household problem given wage and interest rate
    # 3. Obtain labor supply and savings supply
    # 4. Update labor demand = labor supply
    # 5. Update interest rate using bisection
    t_start = time()
    for mm = 1:Opt.MaxIterSteady
        if Opt.display > 0
            println("ITERATION = " * string(mm))
        end

        # I.a Solve worker HJB

        if Opt.display > 0
            println("   Solving worker HJB")
        end
        # solve for labor supply at zero savings

        lab0update!(lab0, TransW, what, Rb, G, P, Opt)

        # iterate on worker HJB until convergence
        iterHJB!(v, Pol, A, B, Rb, what, TransW, lab0, adrift, P, G, Opt, TempW)

        # I.b Solve worker KFE

        if Opt.display > 0
            println("   Solving worker KFE")
        end
        minusATassemble!(minusAT, A, P.zeta, G.iaz, G.jaz)
        linear_solve!(gtemp, minusAT, G.b_vec, TempW,
            Opt.worker_KFE_solver, Opt.bicgstab_tol;
            kfe=true, label="worker KFE")
        g_tilde .= gtemp ./ sum(gtemp)
        g .= reshape(g_tilde, G.N_b, G.N_a, G.N_j, G.N_k)

        # marginal densities
        g_marg_a .= sum(g, dims=[4, 3, 1])[:]

        if P.share_entre > 0.0

            # II.a Solve entrepreneur HJB

            # iterate on entrepreneur HJB until convergence
            if Opt.display > 0
                println("   Solving entrepreneur HJB")
            end
            iterHJBentre!(ve, PolE, Aentre_local, Bentre, RbE, rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants,
                adriftE, P, G, Opt, TempE)

            # II.b Solve entrepreneur KFE

            if Opt.display > 0
                println("   Solving entrepreneur KFE")
            end
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
            ge .= reshape(ge_tilde, G.N_be, G.N_a, G.N_e, G.N_y)

            # marginal densities
            ge_marg_a .= reshape(sum(ge, dims=[4, 3, 1]), G.N_a)
            ge_marg_ey .= reshape(sum(ge, dims=[2, 1]), G.N_e, G.N_y)

        end

        # III. Compute aggregates

        if Opt.display > 0
            println("   Computing aggregates")
        end

        AggrIlliqSavW = dot(g_marg_a, G.a_grid) .* (1.0 .- P.share_entre)

        if P.share_entre > 0
            AggrIlliqSavE = dot(ge_marg_a, G.a_grid) .* P.share_entre
        elseif P.share_entre == 0
            AggrIlliqSavE = 0
        end

        AggrIlliqSav = AggrIlliqSavW .+ AggrIlliqSavE

        for jj = 1:G.N
            LabSupBAJK[jj] = G.z_grid_MAT[jj] .* Pol.lab[jj] .* g[jj]
        end
        AggrLabSup = sum(LabSupBAJK[:]) .* (1.0 .- P.share_entre)
        LabDemE[:] .= LabDemENoLev[:] .* PolE.θ .^ ((P.alphaaE .* P.nu) / (1.0 .- P.nu * (1.0 .- P.alphaaE))) # labor demand of entrepreneurs
        AggrLabDemE = dot(LabDemE[:], ge[:]) .* P.share_entre
        AggrLabDemC = max(AggrLabSup .- AggrLabDemE, 1e-16)
        AggrCapDemE = dot((PolE.θ .- 1.0) .* G.e_grid_MAT[:], ge[:]) .* P.share_entre

        rk = P.alphaa .* mc .* CapLabRatio .^ (P.alphaa .- 1.0)
        AggrCapSup = AggrIlliqSav ./ (1.0 .+ (1.0 .- omega_entre) .* P.omega_liq ./ (rk .- P.delta) .* (1.0 .- mc .- P.theta ./ 2.0 .* P.infl .^ 2.0) ./ AggrCapSup .* Output_total)
        AggrCapDemC = max(AggrCapSup .- AggrCapDemE, 1e-16)
        CapLabRatioNew = AggrCapDemC ./ AggrLabDemC

        ChangeCapLabRatio = CapLabRatioNew ./ CapLabRatio .- 1.0
        if Opt.display > 0
            println("   Average θ = " * string(dot(ge[:], PolE.θ[:])))
            println("   Aggregate leverage (K/E) = " * string(dot(ge[:], PolE.θ[:] .* G.e_grid_MAT[:]) ./ dot(ge[:], G.e_grid_MAT[:])))
            println("   Change K/N = " * string(ChangeCapLabRatio) * "x100%")
            println("   Old K/N = " * string(CapLabRatio))
            println("   New K/N = " * string(CapLabRatioNew))
        end
        if isnan(ChangeCapLabRatio)
            @error "Change in capital-labor-ratio is not a number!"
        end
        CapLabRatio = CapLabRatio + Opt.caplab_stepsize .* (CapLabRatioNew .- CapLabRatio)
        AggrCapDemC = CapLabRatio .* AggrLabDemC
        AggrCapSup = AggrCapDemC .+ AggrCapDemE

        # IV. Update

        if Opt.display > 0
            println("   Updating aggregate variables")
        else
            print(".")
        end
        flush(stdout)

        w = (1.0 .- P.alphaa) .* mc .* CapLabRatio .^ P.alphaa
        what = (rb * P.inflw) * P.thetaw / P.epsilonw + (1.0 + P.tauW) * (1.0 - 1.0 / P.epsilonw) * w  # wage paid to workers
        subs_payment_union = P.tauW * w * AggrLabSup # subsidy paid to union
        Tunion = subs_payment_union  # lump-sum tax on union
        UnionProfits = ((1.0 + P.tauW) * w - what) * AggrLabSup - P.thetaw / 2.0 * P.inflw^2.0 * AggrLabSup - Tunion # Remaining profits of union
        OutputC = AggrLabDemC .^ (1.0 .- P.alphaa) .* AggrCapDemC .^ P.alphaa
        OutputE[:] .= OutputEnoLev[:] .* PolE.θ .^ ((P.alphaaE .* P.nu) ./ (1.0 .- P.nu .* (1.0 .- P.alphaaE))) # output of entrepreneurs
        AggrOutputE = dot(OutputE[:], ge[:]) .* P.share_entre
        Output_total = OutputC .+ AggrOutputE
        omega_entre = AggrOutputE ./ Output_total
        AdjCosts = P.theta ./ 2.0 .* P.infl .^ 2.0 .* Output_total
        rk = P.alphaa .* mc .* CapLabRatio .^ (P.alphaa .- 1.0)
        RealProfits = (1.0 .- pb) .* Output_total .- AdjCosts

        RealProfitsToW = (1.0 .- P.omega_liq) .* RealProfits + UnionProfits # total profits that go to workers liquid account in proportion to labor income
        RealProfitsToE = P.omega_liq .* omega_entre .* RealProfits # total profits that go to entrepreneurs in proportion to output income
        AvgPerCapitaRealProfitsToW = RealProfitsToW ./ (1.0 .- P.share_entre) # average profits PER WORKER going to liquid account
        AvgPerCapitaRealProfitsToE = RealProfitsToE ./ P.share_entre # average profits PER ENTREPRENEURS going to business income

        ge_marg_b .= sum(ge, dims=[4, 3, 2])[:]
        g_marg_b .= sum(g, dims=[4, 3, 2])[:]
        ge_marg_e .= reshape(sum(ge, dims=[4, 2, 1]), G.N_e)
        AggrLiqSav = dot(g_marg_b, G.b_grid) .* (1.0 .- P.share_entre) .+ dot(ge_marg_b, G.b_grid_e) .* P.share_entre
        AggrFirmSav = dot(ge_marg_e, G.e_grid) .* P.share_entre
        Total_assets = AggrLiqSav .+ AggrIlliqSav .+ AggrFirmSav

        zeta_corr = 1.0 .- G.inh_wealth ./ Total_assets # adjust zerta_corr term to account for the share of wealth that is inherited (and not paid out as annuity)

        Rb .= rb .+ P.kappa .* (G.b_grid .< 0) .+ P.zeta .* zeta_corr
        RbE .= rb .+ P.kappa .* (G.b_grid_e .< 0) .+ P.zeta .* zeta_corr
        ra = rk .- P.delta .+ P.zeta .* zeta_corr
        Trans = P.transshare .* Output_total
        TransW .= Trans .+ (1.0 .- P.tau) .* AvgPerCapitaRealProfitsToW .* G.z_grid_MAT ./ G.AvgLabProduct

        if P.share_entre > 0.0
            LabDemENoLev .= (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (1.0 .- P.alphaaE) .* P.nu .* G.e_grid_MAT .^ (P.alphaaE * P.nu) / w) .^ (1.0 ./ (1.0 .- (1.0 .- P.alphaaE) .* P.nu)) # optimal labor demand of entrepreneurs
            OutputEnoLev .= P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (G.e_grid_MAT .^ P.alphaaE .* LabDemENoLev .^ (1.0 .- P.alphaaE)) .^ P.nu
            RealGrossProfitsnoLevPerCap = (pe .* P.entre_prod .* reshape(G.y_grid, 1, 1, 1, G.N_y) .* (P.nu * (1.0 - P.alphaaE) / w)^((1 - P.alphaaE) * P.nu)) .^ (1 / (1 - (1 - P.alphaaE) * P.nu)) * (1 - (1 - P.alphaaE) * P.nu) .* G.e_grid_MAT .^ ((P.nu - 1) / (1 - (1 - P.alphaaE) * P.nu))
            RealGrossProfitsnoLevPerCap[G.e_grid_MAT.==0.0] .= 1.0
            TransE .= Trans .+ (1.0 .- P.tauE) .* AvgPerCapitaRealProfitsToE .* OutputE ./ (AggrOutputE ./ P.share_entre)
        end

        # dampen illiquid returns
        ltau0a = ra .* (G.a_max .* 0.999) .^ (1 - Opt.ltau)
        adriftE = ra .* G.a_grid_MAT_e .- ltau0a .* (G.a_grid_MAT_e .^ Opt.ltau)
        adrift = ra .* G.a_grid_MAT .- ltau0a .* (G.a_grid_MAT .^ Opt.ltau)
        if P.share_entre > 0.0
            FirmIncomeConstants .=(-P.delta_e + P.tauE * expected_depreciation_rate + P.zeta * zeta_corr) .* G.e_grid_MAT # Constant firm income terms: tax deduction for depreciation and annuity (from accidental bequests)
        end

        if abs(ChangeCapLabRatio) < Opt.tol
            timer = string(round(time() .- t_start, digits=1))
            println("Steady state found after " * string(mm) * " iterations and $timer seconds")
            break
        elseif mm == Opt.MaxIterSteady
            timer = string(round(time() .- t_start, digits=1))
            @error "Maximum number of $mm iterations reached after $timer seconds, no steady state found"
            break
        end
    end

    AggrLiqSavW = dot(g_marg_b, G.b_grid) .* (1.0 .- P.share_entre)
    AggrLiqSavE = dot(ge_marg_b, G.b_grid_e) .* P.share_entre
    RealGrossProfitsE[:] .= pe .* OutputE[:] .- w .* LabDemE[:] .- rk .* G.e_grid_MAT[:] .* (PolE.θ .- 1.0)
    TaxableBusinessIncomeE = RealGrossProfitsE .+ AvgPerCapitaRealProfitsToE .* OutputE ./ (AggrOutputE ./ P.share_entre) .- expected_depreciation_rate .* G.e_grid_MAT

    if P.share_entre > 0.0
        opt_firm_size = ((pe * P.entre_prod * G.y_grid * P.nu) * (P.alphaaE / (ra - P.zeta * zeta_corr))^(1 - (1 - P.alphaaE) * P.nu) * ((1 - P.alphaaE) / w)^((1 - P.alphaaE) * P.nu)) .^ (1 / (1 - P.nu))
        ind_opt_firm_size = Vector{Int64}(undef, G.N_y)
        for yy = 1:G.N_y
            ind_opt_firm_size[yy] = argmax(G.e_grid[G.e_grid.<opt_firm_size[yy]]) # optimal firm size for each productivity type
        end
    end

    q = RealProfits * (1.0 - omega_entre) * P.omega_liq ./ (rk .- P.delta)
    BondDem = AggrLiqSav

    Trans = P.transshare .* Output_total # government transfer
    if P.share_entre > 0.0
        tax_rev_entre = P.share_entre .* P.tauE .* dot(TaxableBusinessIncomeE[:], ge[:])
    else
        tax_rev_entre = 0.0
    end
    tax_rev_worker = P.tau .* (w .* AggrLabSup .+ RealProfitsToW)
    subs_payment_union = P.tauW * w * AggrLabSup # subsidy paid to union
    Tunion = subs_payment_union  # lump-sum tax on union
    UnionProfits = ((1.0 + P.tauW) * w - what) * AggrLabSup - P.thetaw / 2.0 * P.inflw^2.0 * AggrLabSup - Tunion # Remaining profits of union
    Gexp = tax_rev_worker .+ tax_rev_entre .- rb .* BondDem .- Trans + Tunion - subs_payment_union # government expenditure

    PolE.f .= TempE.If_B[:] .* TempE.f_B[:] .+ TempE.If_F[:] .* TempE.f_F[:]
    PolE.de .= TempE.Id_B[:] .* TempE.d_B[:] .+ TempE.Id_F[:] .* TempE.d_F[:]
    Pol.d .= TempW.Id_B[:] .* TempW.d_B[:] .+ TempW.Id_F[:] .* TempW.d_F[:]

    InvC_M = P.delta .* AggrCapDemC # market intermediated investment in corporate sector
    InvE_M = P.delta .* AggrCapDemE # market intermediated investment in entrepreneurial sector
    Inv_M = InvC_M + InvE_M # market intermediated investment, should equal P.delta .* AggrCapSup
    OperatingIncomeE = TempE.edrift0 .-
                       P.zeta .* zeta_corr .* G.e_grid_MAT
    OperatingLossE = min.(OperatingIncomeE, 0.0)
    InvE_nonM = P.delta_e .* AggrFirmSav .+
                P.share_entre .* dot(PolE.f .+ OperatingLossE[:], ge[:]) # replacement plus chosen investment, net of operating losses absorbed by business equity
    InvE = InvE_M .+ InvE_nonM # total investment in entrepreneurial sector, market intermediated and non-market intermediated
    Inv = Inv_M .+ InvE_nonM # total investment, market intermediate and non-market intermediated. All investment in the corporate sector is market intermediated.

    if P.share_entre > 0.0
        ConsumpW = dot(g[:], Pol.c[:]) .* (1.0 .- P.share_entre)
        ConsumpE = dot(ge[:], PolE.ce[:]) .* P.share_entre
    elseif P.share_entre == 0.0
        ConsumpW = dot(g[:], Pol.c[:])
        ConsumpE = 0.0
    end
    C = ConsumpW .+ ConsumpE

    borrowing_costs_worker_MAT = -P.kappa .* G.b_grid_MAT .* (G.b_grid_MAT .< 0.0)
    borrowing_costs_entre_MAT = -P.kappa .* G.b_grid_MAT_e .* (G.b_grid_MAT_e .< 0.0)
    adj_costs_a_worker_MAT = P.χ0 .* abs.(Pol.d) .+ P.χ1 .* abs.(Pol.d ./ G.a_grid_MAT_adj[:]) .^ P.χ2 .* G.a_grid_MAT_adj[:]

    adj_costs_entre_e_E = P.χe1 .* abs.(PolE.f ./ G.e_grid_MAT_adj[:]) .^ P.χe2 .* G.e_grid_MAT_adj[:]
    adj_costs_entre_a_E = P.χ0 .* abs.(PolE.de) .+ P.χ1 .* abs.(PolE.de ./ G.a_grid_MAT_e_adj[:]) .^ P.χ2 .* G.a_grid_MAT_e_adj[:]
    PortfAdjW = (1.0 .- P.share_entre) .* dot(g[:], adj_costs_a_worker_MAT[:])
    PortfAdjE = P.share_entre .* dot(ge[:], adj_costs_entre_a_E[:])
    CapAdjE = P.share_entre .* dot(ge[:], adj_costs_entre_e_E[:])
    BorrCostsW = (1.0 .- P.share_entre) .* dot(g[:], borrowing_costs_worker_MAT[:])
    BorrCostsE = P.share_entre .* dot(ge[:], borrowing_costs_entre_MAT[:])
    BorrCosts = BorrCostsW .+ BorrCostsE
    PortfAdj = PortfAdjW .+ PortfAdjE
    KW = AggrCapSup .* AggrIlliqSavW ./ AggrIlliqSav
    KE = AggrCapSup .* AggrIlliqSavE ./ AggrIlliqSav
    InvCE = P.delta * KE
    InvCW = P.delta * KW

    ge_marg_bae .= reshape(sum(ge, dims=[4]), G.N_be, G.N_a, G.N_e)
    ge_marg_ba .= reshape(sum(ge, dims=[4, 3]), G.N_be, G.N_a)
    g_marg_ba .= reshape(sum(g, dims=[4, 3]), G.N_b, G.N_a)

    SSminout = SSminstruct(v, ve, lab0, TransW, TransE, gtemp, gtempE, rb, ra, pb, pe, pp, w, what, Output_total, AggrOutputE, RealProfits,
        AggrCapSup, CapLabRatio, omega_entre, zeta_corr)

    SSout = SSstruct(v, ve, lab0, TransW, TransE, Trans, gtemp, gtempE, rb, ra, pb, pe, pp, w, what, Output_total, AggrOutputE, RealProfits,
        AggrCapSup, KW, KE, CapLabRatio, omega_entre, zeta_corr, LabDemE, AggrCapDemE, AggrCapDemC, OutputE,
        RealGrossProfitsE, TaxableBusinessIncomeE, g, ge, ge_marg_b, ge_marg_a, ge_marg_e, ge_marg_ey, ge_marg_ba,
        ge_marg_bae, g_marg_b, g_marg_a, g_marg_ba, AggrLiqSavW, AggrIlliqSavW,
        AggrLiqSavE, AggrIlliqSavE, AggrFirmSav, AggrLabDemE, AggrLabSup, AggrLiqSav,
        AggrIlliqSav, AggrLabDemC, rk, OutputC, AdjCosts, PortfAdjW,
        PortfAdjE, CapAdjE, BorrCostsW, BorrCostsE, BorrCosts, PortfAdj,
        BondDem, q, tax_rev_entre, tax_rev_worker, subs_payment_union, Gexp, Inv, Inv_M, InvC_M, InvE_M, InvE_nonM, InvCW, InvCE, InvE, ConsumpE,
        ConsumpW, C, RealProfitsToW, RealProfitsToE, AvgPerCapitaRealProfitsToW, AvgPerCapitaRealProfitsToE, UnionProfits, Tunion)

    # Some checks on A, Anetre, g, ge
    Aentrerowssumsextrema =extrema(vec(sum(Aentre, dims=2)))
    if Aentrerowssumsextrema[1] < -1e-8 || Aentrerowssumsextrema[2] > 1e-8
        @warn "Row sums of Aentre are not close to zero, extrema are: " * string(Aentrerowssumsextrema)
    end
    Arowsumsextrema = extrema(vec(sum(A, dims=2)))
    if Arowsumsextrema[1] < -1e-8 || Arowsumsextrema[2] > 1e-8
        @warn "Row sums of A are not close to zero, extrema are: " * string(Arowsumsextrema)
    end
    Adiagmax = maximum(vec(diag(A)))
    if Adiagmax > 1e-8
        @warn "Diagonal elements of A are not all negative, maximum is: " * string(Adiagmax)
    end
    Aentrediagmax = maximum(vec(diag(Aentre)))
    if Aentrediagmax > 1e-8
        @warn "Diagonal elements of Aentre are not all negative, maximum is: " * string(Aentrediagmax)
    end
    ia, ja, va = findnz(Aentre)
    Aentreoffdiagmax = minimum(va[ia .!= ja])
    if Aentreoffdiagmax < -1e-8
        @warn "Off-diagonal elements of Aentre are not all non-negative, minimum is: " * string(Aentreoffdiagmax)
    end
    ia, ja, va = findnz(A)
    Aoffdiagmax = minimum(va[ia .!= ja])
    if Aoffdiagmax < -1e-8
        @warn "Off-diagonal elements of A are not all non-negative, minimum is: " * string(Aoffdiagmax)
    end

    gemin = minimum(ge)
    if gemin < -1e-8
        @warn "Some elements of ge are negative, minimum is: " * string(gemin)
        @warn "Clipping ge to be non-negative"
        ge[ge .< 0.0] .= 0.0
        ge .= ge ./ sum(ge) # renormalize
    end
    gmin = minimum(g)
    if gmin < -1e-8
        @warn "Some elements of g are negative, minimum is: " * string(gmin)
        @warn "Clipping g to be non-negative"
        g[g .< 0.0] .= 0.0
        g .= g ./ sum(g) # renormalize
    end

    return SSminout, SSout, Pol, PolE, A, Aentre

end

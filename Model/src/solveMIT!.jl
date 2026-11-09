"`solveMIT!(MIT, MaxAbsCapDev, MaxAbsBondDev, TransPolW, TransPolE, gstart, gestart, SSallend, mp_shock_t, TempW, TempE, A0fix, G, P, Opt, ModelName, doplot, dosave, entre_exp, worker_exp,TransPolWNA=TransPolW,gstartNA=false,shareNA=0.0)`
Solves MIT shock from initial distributions `gstart` and `gestart` with final steady state `SSallend`.
Modifies `MIT` and updates convergence criteria (`MaxAbsCapDev`,`MaxAbsBondDev`) .
`A0fix` = true to keep initial distributions unaffected by stock price q.
See also solveMIT for the version that initializes MIT and convergence criteria"
function solveMIT!(MIT, MaxAbsCapDev, MaxAbsBondDev, TransPolW, TransPolE, gstart, gestart, SSallend, mp_shock_t, TempW, TempE, A0fix, G, P, Opt, ModelName, doplot, dosave, entre_exp, worker_exp, TransPolWNA=TransPolW, gstartNA=false, shareNA=0.0; n_generators::Integer=0, validation_signature=nothing)
    println("Computing MIT transition for " * ModelName)
    if entre_exp != 0
        @warn "Solving with alternative assumption on entrepreneurial behavior (entre_exp = " * string(entre_exp) * ")! Set entre_exp = 0 for true model."
    end
    if worker_exp != 0
        @warn "Solving with alternative assumption on worker behavior (worker_exp = " * string(worker_exp) * ")! Set worker_exp = 0 for true model."
    end

    # Initialize
    MIT.mp_shock_t .= mp_shock_t
    MIT.q_t[end] = SSallend.q
    MIT.mc_t[end] = SSallend.pe
    AggrOutputE_new = Array{Float64}(undef, G.N_t)  # aggregate output of entrepreneurs
    AggrLabDemE_new = Array{Float64}(undef, G.N_t)  # aggregate labor demand of entrepreneurs
    tax_rev_entre_new = Array{Float64}(undef, G.N_t)
    Capital_new = Array{Float64}(undef, G.N_t)
    rb_new = Array{Float64}(undef, G.N_t)
    MIT.wagemarkdown_t .= zeros(G.N_t)
    MIT.what_t .= MIT.wagemarkdown_t .+ SSallend.what
    g = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    g_tilde = Array{Float64}(undef, G.N)
    ge = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    g_tilde_e = Array{Float64}(undef, G.Ne)
    ge_marg_a = Array{Float64}(undef, G.N_a)
    g_marg_a = Array{Float64}(undef, G.N_a)

    if worker_exp == 4 # some workers do no respond:
        g_tildeNA = Array{Float64}(undef, G.N)
        ANA_t = sparse_value_path(G.Az, G.N_t)

        g_marg_aNA = Array{Float64}(undef, G.N_a)
        gNA = deepcopy(gstartNA) #Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    end

    Ainit = deepcopy(G.Aexo)
    # The generators have fixed CSC structure; only their nzval arrays differ by date.
    A_t = sparse_value_path(G.Az, G.N_t)
    Aentre_local_t = sparse_value_path(Ainit, G.N_t)

    if doplot
        nplots = 30
        plots_short = Vector{Plots.Plot{Plots.GRBackend}}(undef, nplots)
        plots_long = Vector{Plots.Plot{Plots.GRBackend}}(undef, nplots)
        for nn = 1:nplots
            plots_long[nn] = plot()
            plots_short[nn] = plot()
        end
    end

    ge_marg_e = vec(sum(gestart, dims=(1, 2, 4)))
    MIT.AggrFirmSav_t[1] = dot(G.e_grid, ge_marg_e) .* P.share_entre

    for it = (MIT.it[1]+1):Opt.MaxIterMIT
        t_start = time()
        MIT.OutputC_t .= MIT.AggrCapDemC_t .^ P.alphaa .* MIT.AggrLabDemC_t .^ (1.0 .- P.alphaa)
        MIT.Output_total_t .= MIT.OutputC_t .+ MIT.AggrOutputE_t
        MIT.omega_entre .= MIT.AggrOutputE_t ./ MIT.Output_total_t
        MIT.dinfl_t .= diff(MIT.infl_t)
        MIT.dOutput_t .= diff(MIT.Output_total_t)
        MIT.dinflw_t .= diff(MIT.inflw_t)
        MIT.dLab_t .= diff(MIT.AggrLabSup_t)
        mdadj = 0.5
        MIT.wagemarkdown_t[end] = (MIT.rb_t[end] .* MIT.inflw_t[end]) .* P.thetaw ./ P.epsilonw * mdadj .+ (1.0 - mdadj) .* MIT.wagemarkdown_t[end]
        MIT.wagemarkdown_t[1:(end-1)] .= (((MIT.rb_t[1:(end-1)] .- MIT.dLab_t ./ G.dt[1:(end-1)] ./ MIT.AggrLabSup_t[1:(end-1)]) .* MIT.inflw_t[1:(end-1)] .- MIT.dinflw_t ./ G.dt[1:(end-1)]) .* P.thetaw ./ P.epsilonw) .* mdadj .+ (1.0 - mdadj) * MIT.wagemarkdown_t[1:(end-1)] # wage paid to workers
        MIT.mc_t[1:(end-1)] .= SSallend.pe .+ P.theta ./ P.epsilon .* (((MIT.rk_t[1:(end-1)] .- P.delta) .- (MIT.dOutput_t ./ G.dt[1:(end-1)]) ./ MIT.Output_total_t[1:(end-1)]) .* MIT.infl_t[1:(end-1)] .- MIT.dinfl_t ./ G.dt[1:(end-1)])

        MIT.pe_t .= MIT.mc_t
        MIT.rk_t .= P.alphaa .* MIT.mc_t .* (MIT.AggrCapDemC_t ./ MIT.AggrLabDemC_t) .^ (P.alphaa .- 1.0)
        MIT.w_t .= (1.0 .- P.alphaa) .* MIT.mc_t .* (MIT.AggrCapDemC_t ./ MIT.AggrLabDemC_t) .^ P.alphaa
        whatadj = 1.0
        MIT.what_t .= ((1.0 .+ P.tauW) .* (1.0 .- 1.0 ./ P.epsilonw) .* MIT.w_t .+ MIT.wagemarkdown_t) .* whatadj .+ (1.0 .- whatadj) .* MIT.what_t
        MIT.UnionProfits_t .= ((1.0 .+ P.tauW) .* MIT.w_t .- MIT.what_t) .* MIT.AggrLabSup_t .- P.thetaw ./ 2.0 .* MIT.inflw_t .^ 2.0 .* MIT.AggrLabSup_t .- MIT.Tunion_t # Remaining profits of union
        MIT.AdjCosts_t .= P.theta ./ 2.0 .* MIT.infl_t .^ 2.0 .* MIT.Output_total_t
        MIT.RealProfits_t .= (1.0 .- MIT.mc_t) .* MIT.Output_total_t .- MIT.AdjCosts_t
        for tt = (G.N_t-1):-1:1
            MIT.q_t[tt] = (MIT.q_t[tt+1] .+ MIT.RealProfits_t[tt] .* (1.0 .- MIT.omega_entre[tt]) .* P.omega_liq .* G.dt[tt]) ./ (1.0 .+ G.dt[tt] .* (MIT.rk_t[tt] .- P.delta))

        end
        MIT.ra_t .= MIT.rk_t .- P.delta .+ P.zeta .* MIT.zeta_corr_t

        # Government lump-sum transfer with constant spending; Transfers adjust!
        MIT.subs_payment_union_t .= P.tauW .* MIT.w_t .* MIT.AggrLabSup_t
        MIT.Trans_t .= P.tau .* (MIT.w_t .* MIT.AggrLabSup_t .+ MIT.RealProfitsToW_t) .- MIT.rb_t .* SSallend.BondDem .+ MIT.tax_rev_entre_t .- MIT.G_t .+ MIT.Tunion_t .- MIT.subs_payment_union_t

        # Compute aggregate illiquid savings from distributions gstart and gestart
        ge_marg_a .= vec(sum(gestart, dims=(1, 3, 4)))
        g_marg_a .= vec(sum(gstart, dims=(1, 3, 4)))
        AggrIlliqSavW = dot(g_marg_a, G.a_grid) .* (1.0 .- P.share_entre)
        AggrIlliqSavE = dot(ge_marg_a, G.a_grid) .* P.share_entre
        AggrIlliqSav = AggrIlliqSavW .+ AggrIlliqSavE
        if worker_exp == 4
            g_marg_aNA .= vec(sum(gstartNA, dims=(1, 3, 4)))
            AggrIlliqSavWNA = dot(g_marg_aNA, G.a_grid) .* (1.0 .- P.share_entre)
            AggrIlliqSavW = AggrIlliqSavW * (1.0 - shareNA) + AggrIlliqSavWNA * shareNA
            AggrIlliqSav = AggrIlliqSavE + AggrIlliqSavW
        end
        # Compute drop in A due to change in q and adjust distribution accordingly or not, depending on whether we want to enforce a given distribution (convergence to steady state without shock) or not (shock)
        if A0fix # do not change initial distribution depending on q
            ge .= gestart
            g .= gstart
            g_tilde_e .= vec(gestart)
            g_tilde .= vec(gstart)
            if worker_exp == 4
                gNA .= gstartNA
                g_tildeNA .= vec(gNA)
            end
        else # adjust distribution depending on q
            # 1. Compute how much the illiquid asset position of each household changes in period 1 because of period 1 change of q: A_0_drop
            if mp_shock_t[1] <= 0.0 # expansionary shock
                A_0_drop = min((MIT.q_t[1] .+ MIT.AggrCapSup_t[1]) ./ AggrIlliqSav, 1.0)
            elseif mp_shock_t[1] > 0.0 # contractionary shock
                A_0_drop = max((MIT.q_t[1] .+ MIT.AggrCapSup_t[1]) ./ AggrIlliqSav, 1.0)
            end
            # 2. Adjust initial distributions accordingly
            g .= gstart
            if A_0_drop != 1.0
                initial_q_change!(g, A_0_drop, G.a_grid, G.N_a, 2)  # adjust the initial density for proportional drop in q
            end
            g_tilde .= vec(g)
            if worker_exp == 4
                gNA .= gstartNA
                if A_0_drop != 1.0
                    initial_q_change!(gNA, A_0_drop, G.a_grid, G.N_a, 2)  # adjust the initial density for proportional drop in q
                end
                g_tildeNA .= vec(gNA)
            end
            if P.share_entre > 0.0
                ge .= gestart
                if A_0_drop != 1.0
                    initial_q_change!(ge, A_0_drop, G.a_grid, G.N_a, 2)
                end
                g_tilde_e .= vec(ge)
            else
                ge .= 0.0
                g_tilde_e .= 0.0
            end
        end

        ge_marg_a .= vec(sum(ge, dims=(1, 3, 4)))
        g_marg_a .= vec(sum(g, dims=(1, 3, 4)))
        MIT.AggrIlliqSavW_t[1] = dot(g_marg_a, G.a_grid) .* (1.0 .- P.share_entre)
        MIT.AggrIlliqSavE_t[1] = dot(ge_marg_a, G.a_grid) .* P.share_entre
        MIT.AggrIlliqSav_t[1] = MIT.AggrIlliqSavW_t[1] .+ MIT.AggrIlliqSavE_t[1]
        if worker_exp == 4
            g_marg_aNA .= vec(sum(gNA, dims=(1, 3, 4)))
            AggrIlliqSavWNA = dot(g_marg_aNA, G.a_grid) .* (1.0 .- P.share_entre)
            MIT.AggrIlliqSavW_t[1] = MIT.AggrIlliqSavW_t[1] * (1.0 - shareNA) + AggrIlliqSavWNA * shareNA
            MIT.AggrIlliqSav_t[1] = MIT.AggrIlliqSavW_t[1] + MIT.AggrIlliqSavE_t[1]
        end

        if worker_exp == 4
            MITupdate!(MIT, AggrOutputE_new, AggrLabDemE_new, tax_rev_entre_new, g_tilde, g_tilde_e, SSallend,
                A_t, Aentre_local_t, TransPolW, TransPolE, g, ge, TempW, TempE, P, G, Opt, entre_exp, worker_exp, g_tildeNA, ANA_t, TransPolWNA, gNA, shareNA)
        else
            MITupdate!(MIT, AggrOutputE_new, AggrLabDemE_new, tax_rev_entre_new, g_tilde, g_tilde_e, SSallend,
                A_t, Aentre_local_t, TransPolW, TransPolE, g, ge, TempW, TempE, P, G, Opt, entre_exp, worker_exp)
        end

        if P.share_entre > 0
            # update aggregate variables from entrepr. sector conservatively
            MIT.AggrOutputE_t .= Opt.Entreadjust_MIT .* AggrOutputE_new .+ (1.0 .- Opt.Entreadjust_MIT) .* MIT.AggrOutputE_t
            MIT.AggrLabDemE_t .= Opt.Entreadjust_MIT .* AggrLabDemE_new .+ (1.0 .- Opt.Entreadjust_MIT) .* MIT.AggrLabDemE_t
            MIT.tax_rev_entre_t .= Opt.Entreadjust_MIT .* tax_rev_entre_new .+ (1.0 .- Opt.Entreadjust_MIT) .* MIT.tax_rev_entre_t
        end

        MIT.AggrLabDemC_t .= MIT.AggrLabSup_t .- MIT.AggrLabDemE_t
        MIT.AggrCapDemC_t .= MIT.AggrCapSup_t .- MIT.AggrCapDemE_t
        MIT.AggrIlliqSav_t[2:end] .= MIT.AggrIlliqSavW_t[2:end] .+ MIT.AggrIlliqSavE_t[2:end]
        MIT.AggrLiqSav_t[2:end] .= MIT.AggrLiqSavW_t[2:end] .+ MIT.AggrLiqSavE_t[2:end]

        Capital_new .= MIT.AggrIlliqSav_t[1:(end-1)] .- MIT.q_t[1:end]
        MaxAbsCapDev[it] = maximum(abs.(MIT.AggrCapSup_t .- Capital_new) ./ SSallend.Output_total)
        MaxAbsBondDev[it] = maximum(abs.(MIT.AggrLiqSav_t .- SSallend.BondDem) ./ SSallend.Output_total)

        MIT.AggrCapSup_t .= Opt.stepstickytransK .* Capital_new .+ (1.0 .- Opt.stepstickytransK) .* MIT.AggrCapSup_t
        MIT.AggrCapDemC_t .= MIT.AggrCapSup_t .- MIT.AggrCapDemE_t
        MIT.KW_t .= MIT.AggrCapSup_t .* MIT.AggrIlliqSavW_t[1:(end-1)] ./ MIT.AggrIlliqSav_t[1:(end-1)]
        MIT.KE_t .= MIT.AggrCapSup_t .* MIT.AggrIlliqSavE_t[1:(end-1)] ./ MIT.AggrIlliqSav_t[1:(end-1)]
        for tt = (G.N_t-1):-1:1
            MIT.InvC_M_t[tt] = (MIT.AggrCapDemC_t[tt+1] .- MIT.AggrCapDemC_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT.AggrCapDemC_t[tt]
            MIT.InvE_M_t[tt] = (MIT.AggrCapDemE_t[tt+1] .- MIT.AggrCapDemE_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT.AggrCapDemE_t[tt]
            MIT.Inv_M_t[tt] = MIT.InvC_M_t[tt] .+ MIT.InvE_M_t[tt]
            MIT.InvCW_t[tt] = (MIT.KW_t[tt+1] .- MIT.KW_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT.KW_t[tt]
            MIT.InvCE_t[tt] = (MIT.KE_t[tt+1] .- MIT.KE_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT.KE_t[tt]
        end
        MIT.InvE_t .= MIT.InvE_M_t .+ MIT.InvE_nonM_t
        MIT.Inv_t .= MIT.Inv_M_t .+ MIT.InvE_nonM_t
        MIT.zeta_corr_t .= 1.0 .- G.inh_wealth ./ (MIT.AggrLiqSav_t[1:(end-1)] .+ MIT.AggrIlliqSav_t[1:(end-1)] .+ MIT.AggrFirmSav_t[1:(end-1)]) # adjust zerta_corr term to account for the share of wealth that is inherited (and not paid out as annuity)

        MIT.Excess_Sup_Liq_t .= MIT.AggrLiqSav_t .- SSallend.BondDem
        rb_new .= MIT.rb_t .- (diff(MIT.Excess_Sup_Liq_t) ./ (Opt.bondadjust_MIT .* G.dt) .+ MIT.Excess_Sup_Liq_t[1:(end-1)]) ./ SSallend.Output_total
        MIT.rb_t .= Opt.stepstickytransB .* rb_new .+ (1.0 .- Opt.stepstickytransB) .* MIT.rb_t

        MIT.infl_t .= (MIT.rb_t .- mp_shock_t .- SSallend.rb) ./ (P.phi .- 1.0)
        MIT.inflw_t[1:(end-1)] .= diff(MIT.w_t) ./ G.dt[1:(end-1)] ./ MIT.w_t[1:(end-1)] .+ MIT.infl_t[1:(end-1)]

        MIT.rb_t[end] = SSallend.rb
        MIT.infl_t[end] = (MIT.rb_t[end] - mp_shock_t[end] - SSallend.rb) / (P.phi - 1.0)

        if Opt.display > 0
            println("ITERATION = " * string(it) * ", time for iteration = " * string(round(time() .- t_start, digits=1)) * " sec., Max |Excess Illiquid Savings| = " * string(MaxAbsCapDev[it]) * ", Max |Excess Liquid Savings| = " * string(MaxAbsBondDev[it]))
        end

        # plot some intermediate outcomes
        if doplot
            plotMITconv(MaxAbsCapDev, MaxAbsBondDev, G.dt_cum, Capital_new, MIT.AggrCapSup_t, it, ModelName; output_dir=joinpath(Opt.outpath, "Diagnostics"))

            plotMIT!(plots_long, MIT, mp_shock_t, SSallend, 999, G, ModelName; output_dir=joinpath(Opt.outpath, "Diagnostics"))
            plotMIT!(plots_short, MIT, mp_shock_t, SSallend, 5, G, ModelName; output_dir=joinpath(Opt.outpath, "Diagnostics"))

        end
        MIT.it[1] = it
        if dosave
            @save Opt.outpath * "MIT_policies_" * string(ModelName) * ".jld2" TransPolW TransPolE validation_signature #A_t Aentre_t
            @save Opt.outpath * "MIT_" * string(ModelName) * ".jld2" MIT SSallend P Opt MaxAbsCapDev MaxAbsBondDev validation_signature
        end

        if (MaxAbsCapDev[it] < Opt.MIT_tol) & (MaxAbsBondDev[it] < Opt.MIT_tol)
            println("Converged!")
            MIT.converged[1] = true
            MIT.crit .= [MaxAbsCapDev[it], MaxAbsBondDev[it]]
            break
        else
            MIT.converged[1] = false
        end
    end

    if dosave
        @save Opt.outpath * "MIT_policies_" * string(ModelName) * ".jld2" TransPolW TransPolE validation_signature #A_t Aentre_t
        @save Opt.outpath * "MIT_" * string(ModelName) * ".jld2" MIT SSallend P Opt MaxAbsCapDev MaxAbsBondDev validation_signature
    end

    println("Done")

    if n_generators == 0
        return nothing
    end
    # Construct full generator including leverage-dependent risky jumps for entrepreneurs
    entrepreneur_generators = [complete_entre_generator(Aentre_local_t[n], TransPolE[n].θ, P, G, Opt) for n = 1:n_generators]
    return (workers=A_t[1:n_generators],entrepreneurs=entrepreneur_generators,)
end

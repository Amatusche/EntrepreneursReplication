function main_fct()

    Opt = options()
    println("Leverage = ", Opt.leverage_mode)
    P = parameters(Opt)
    G = grids(P, Opt)

    if Opt.stickywages
        @set! P.theta = 0.0
        @set! P.epsilon = Inf
        @set! P.epsilonw = 100.0
        @set! P.thetaw = 10.0
        @set! P.tauW = 1.0 ./ (P.epsilonw .- 1.0)
    end

    mkpath(Opt.outpath)
    SSInitGuess = steady_init(G, P, Opt)
    try SSInitGuess = load(Opt.outpath * "SS.jld2", "SS") catch end

    SS, SSall, Pol, PolE, A, Aentre = findsteady(SSInitGuess, P, G, Opt)
    SSStats = stats(SSall, PolE, Aentre, G, P, Opt)
    @save Opt.outpath * "SS.jld2" SS SSall SSStats Pol PolE Aentre P G Opt

    if Opt.calibrate
        P, Targets, Stats = calibrate(SSInitGuess, P, Opt)
        G = grids(P, Opt)
        SS, SSall, Pol, PolE, A, Aentre = findsteady(SSInitGuess, P, G, Opt)
        SSStats = stats(SSall, PolE, Aentre, G, P, Opt)
        @save Opt.outpath * "SS.jld2" SS SSall SSStats Pol PolE Aentre P G Opt
    end

    MITInitGuess, TransPolWInitGuess, TransPolEInitGuess = MIT_init(SSall, G.mp_shock_t, Pol, PolE, G, P)
    try MITInitGuess = load(Opt.outpath * "MIT_Base.jld2", "MIT") catch end
    try TransPolWInitGuess, TransPolEInitGuess = load(Opt.outpath * "MIT_policies_Base.jld2", "TransPolW", "TransPolE") catch end

    doplot = true
    dosave = true
    entre_exp = 0 # no counterfactual behavior of entrepreneurs
    worker_exp = 0
    MIT, TransPolW, TransPolE, MaxAbsCapDev, MaxAbsBondDev = solveMIT(MITInitGuess, TransPolWInitGuess, TransPolEInitGuess, SSall.g, SSall.ge, SSall, G.mp_shock_t, false, G, P, Opt, "Base", doplot, dosave, entre_exp, worker_exp, copy_inputs=false)
    Partial = partial(MIT, TransPolW, TransPolE, SSall, ["rb"], [MIT.rb_t], P, G, Opt)
    MITStats = statsMIT(MIT, TransPolW, TransPolE, nothing, Partial, nothing, SSall, Pol, PolE, A, Aentre, G, P, Opt)
    @save Opt.outpath * "MITStats_Base.jld2" MITStats Partial

    entre_exp = 3 # only rb perceived unchanged for entrepreneurs' decisions, drift according to actual prices, i.e., they are surprised
    MITInitGuess, TransPolWInitGuess, TransPolEInitGuess = MIT_init(SSall, G.mp_shock_t, Pol, PolE, G, P)
    try MITInitGuess = load(Opt.outpath * "MIT_Base.jld2", "MIT") catch end
    try TransPolWInitGuess, TransPolEInitGuess = load(Opt.outpath * "MIT_policies_Base.jld2", "TransPolW", "TransPolE") catch end
    MITInitGuess.it[1] = 0
    MITInitGuess.converged[1] = false
    MITInitGuess.crit .= Inf
    try MITInitGuess = load(Opt.outpath*"MIT_entre_exp3.jld2", "MIT") catch end
    try TransPolWInitGuess, TransPolEInitGuess = load(Opt.outpath * "MIT_policies_entre_exp3.jld2", "TransPolW", "TransPolE") catch end
    MIT_entre_exp3, TransPolW_entreexp3, TransPolE_entreexp3, MaxAbsCapDev, MaxAbsBondDev = solveMIT(MITInitGuess, TransPolWInitGuess, TransPolEInitGuess, SSall.g, SSall.ge, SSall, G.mp_shock_t, false, G, P, Opt, "entre_exp3", doplot, dosave, entre_exp, worker_exp)

    worker_exp = 3
    entre_exp = 0
    MITInitGuess, TransPolWInitGuess, TransPolEInitGuess = MIT_init(SSall, G.mp_shock_t, Pol, PolE, G, P)
    try MITInitGuess = load(Opt.outpath * "MIT_Base.jld2", "MIT") catch end
    try TransPolWInitGuess, TransPolEInitGuess = load(Opt.outpath * "MIT_policies_Base.jld2", "TransPolW", "TransPolE") catch end
    MITInitGuess.it[1] = 0
    MITInitGuess.converged[1] = false
    MITInitGuess.crit .= Inf
    try MITInitGuess = load(Opt.outpath*"MIT_worker_exp3.jld2", "MIT") catch end
    try TransPolWInitGuess, TransPolEInitGuess = load(Opt.outpath * "MIT_policies_worker_exp3.jld2", "TransPolW", "TransPolE") catch end
    MIT_worke_exp3, TransPolW_workerexp3, TransPolE_workerexp3, MaxAbsCapDev, MaxAbsBondDev = solveMIT(MITInitGuess, TransPolWInitGuess, TransPolEInitGuess, SSall.g, SSall.ge, SSall, G.mp_shock_t, false, G, P, Opt, "worker_exp3", doplot, dosave, entre_exp, worker_exp)

    Optbeq = deepcopy(Opt)
    beqsize_vec = [500.0, 2000.0, 10000.0] # in 1000 USD
    beqfrac_vec = [100.0, 25.0, 5.0] ./ 100.0
    SSbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    PolEbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    Polbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    Aentrebeq_vec = Vector{Any}(undef, length(beqsize_vec))
    MITbeq_vec = Vector{MITstruct}(undef, length(beqsize_vec))
    TransPolWbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    TransPolEbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    SSStatsbeq_vec = Vector{SteadyStats}(undef, length(beqsize_vec))
    MITStatsbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    Partialbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    Pbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    Gbeq_vec = Vector{Any}(undef, length(beqsize_vec))
    name_vec = Vector{String}(undef, length(beqsize_vec))
    for ii = eachindex(beqsize_vec)
        Pbeq = deepcopy(P)
        target_bequest_model_units = beqsize_vec[ii] * 1e3 / Pbeq.NomGDPFac
        @set! Pbeq.inheritance_index_e =
            fractional_grid_index(G.e_grid, target_bequest_model_units)
        @set! Pbeq.frac_inh_wealth = beqfrac_vec[ii] # share of entrepreneurs receiving inheritance
        Pbeq_vec[ii] = Pbeq
        name_vec[ii] = "Beq" * string(beqsize_vec[ii]) * string(Pbeq_vec[ii].frac_inh_wealth)
        Gbeq_vec[ii] = grids(Pbeq_vec[ii], Optbeq)
        implied_bequest_model_units = Gbeq_vec[ii].inh_wealth /
            (Pbeq.share_entre * Pbeq.frac_inh_wealth)
        @assert isapprox(implied_bequest_model_units, target_bequest_model_units;
            rtol=1e-12, atol=1e-12) "Grid-implied bequest does not match requested bequest for regime $(name_vec[ii])"
    end

    beq_results = pmap(
        (Gbeq, Pbeq, beqsize, name) ->
            bequestregimes!(SSall, Pbeq, Gbeq, Optbeq, beqsize, name,
                SSall.AggrLiqSav),
        Gbeq_vec, Pbeq_vec, beqsize_vec, name_vec)

    for ii = eachindex(beqsize_vec)
        result = beq_results[ii]
        SSbeq_vec[ii] = result.SSallbeq
        Polbeq_vec[ii] = result.Polbeq
        PolEbeq_vec[ii] = result.PolEbeq
        Aentrebeq_vec[ii] = result.Aentrebeq
        MITbeq_vec[ii] = result.MITbeq
        TransPolWbeq_vec[ii] = result.TransPolWbeq
        TransPolEbeq_vec[ii] = result.TransPolEbeq

        Partialbeq_vec[ii] =
            partial(MITbeq_vec[ii], TransPolWbeq_vec[ii], TransPolEbeq_vec[ii],
                SSbeq_vec[ii], ["rb"],
                [MIT.rb_t .- SSall.rb .+ SSbeq_vec[ii].rb], Pbeq_vec[ii],
                Gbeq_vec[ii], Optbeq)
        MITStatsbeq_vec[ii] =
            statsMIT(MITbeq_vec[ii], TransPolWbeq_vec[ii], TransPolEbeq_vec[ii], nothing,
                Partialbeq_vec[ii], nothing, SSbeq_vec[ii], Polbeq_vec[ii],
                PolEbeq_vec[ii], beq_results[ii].Abeq, Aentrebeq_vec[ii],
                Gbeq_vec[ii], Pbeq_vec[ii], Optbeq)
        SSStatsbeq_vec[ii] = stats(SSbeq_vec[ii], PolEbeq_vec[ii], Aentrebeq_vec[ii], Gbeq_vec[ii], Pbeq_vec[ii], Optbeq)
    end

    beq_stats1, beq_stats2 = beqStats(SSall, SSStats, MIT, MITStats, Partial, SSbeq_vec, SSStatsbeq_vec, MITbeq_vec, MITStatsbeq_vec, Partialbeq_vec, beqsize_vec, beqfrac_vec)
    @save Opt.outpath*"bequestregimes.jld2" SSbeq_vec MITbeq_vec SSStatsbeq_vec Partialbeq_vec MITStatsbeq_vec beqsize_vec beqfrac_vec beq_stats1 beq_stats2

    shocks = collect(-0.02:0.005:0.02)
    println("###############")
    println("###############")
    println("Computing Phillips curve for baseline!")
    println("###############")
    println("###############")
    MITVec, MaxAbsCapDevVec, MaxAbsBondDevVec =
        solve_phillips_curve_grid(SSall, Pol, PolE, G, P, Opt;
            shocks=shocks, experiment="baseline", output_path=Opt.outpath * "PC.jld2")

    jj = 3
    println("###############")
    println("###############")
    println("Computing Phillips curve at bequest regime " * string(jj) * "!")
    println("###############")
    println("###############")
    path = Opt.outpath * "PC_beq" * string(jj) * ".jld2"
    MITbeq_PC, MaxAbsCapDev_vec, MaxAbsBondDev_vec =
        solve_phillips_curve_grid(
            SSbeq_vec[jj], Polbeq_vec[jj], PolEbeq_vec[jj], Gbeq_vec[jj],
            Pbeq_vec[jj], Optbeq; shocks=shocks,
            name_prefix="beqreg" * string(jj) * "_",
            experiment="bequest_regime_" * string(jj), output_path=path)

end

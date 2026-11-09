function calibrate(SS, P, Opt)
    println("Doing calibration")
    if Opt.calibonlyE
        @warn "Calibrating only entrepreneur problem (no GE)!"
    end

    # Set up targets for internal calibration

    Targets = calibration_targets()

    # Initial guess for parameters. Calibrate both the stationary share and the
    # transition speed of the statistical entrepreneurial-income states.
    estvec = [true; true; true; true; false; true; true; true; true; false; false; false; false; false; true; true]
    x0 = log.([P.ρ*100; P.ρE/P.ρ; P.kappa*100; P.χ1; P.χ2; P.entre_prod; P.log_y_spread; P.sigma_e*10; P.χe1; P.χe2; P.alphaaE; (P.delta_e+(1-P.eta/(P.eta+1))*P.lambda_e*Opt.entreriskjump)*10; exp.(P.lambda_e); exp.(P.eta); P.avgtransrate; P.lowstateshare])[estvec];

    #           ρ           ρE        κ             χ1        χ2  entre_prod   log_y_spread    σe          χe1    χe2  avgtransrate

    weights = [true, true, false, true, false, true, false, false, true, true, false, false, false, true, false, true, true, true]
    #          A/Y  B/Y  Phtm Whtm negLiqShare Lab at e Liq E Illiq E Htm E NWratio K/Y Htmshare Emp_ratio Immobility Rank autocorr. Large-loss prob. Large-gain prob. Avg. large-loss size

    upper = log.([0.021*100; 1.50; 0.040*100; 1.25; 2.0; 2.5; 1.60; 0.160*10; 1.50; 2.0; 0.4; (P.delta_e+(1-P.eta/(P.eta+1))*P.lambda_e*Opt.entreriskjump).*20.0; exp(0.050); exp(10.0); 0.20; 0.30])[estvec];
    lower = log.([0.0155*100; 0.95; 0.015*100; 0.25; 1.0; 1.85; 0.50; 0.075*10; 0.05; 1.0; 0.2; (P.delta_e+(1-P.eta/(P.eta+1))*P.lambda_e*Opt.entreriskjump)*10.0; exp(0.005); exp(3.5); 0.040; 0.03])[estvec];
    logs = [true; true; true; true; true; true; true; true; true; true; true; true; false; false; true; true][estvec];
    expon = [2; 0; 2; 0; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0; 0][estvec];

    ParamNames = ["ρ", "ρE/ρ", "κ", "χ1", "χ2", "entre_prod", "log_y_spread", "σe", "χe1", "χe2", "alphaaE", "delta_e", "lambda_e", "eta", "avgtransrate", "lowstateshare"][estvec]
    TargetsNames = ["A/Y", "B/Y", "Phtm", "Whtm", "negLiqShare", "Lab at e", "Liq E", "Illiq E", "Htm E", "NWratio", "K/Y", "Htmshare", "Emp_Necessity/Emp_Opportunity", "EInc_immobility_ratio", "EInc_rank_autocorrelation", "EInc_large_loss_probability", "EInc_large_gain_probability", "EInc_avg_deciles_lost"][weights]
    println(string(length(ParamNames) ) * " parameters calibrated: ", ParamNames)
    println(string(length(TargetsNames)) * " targets: ", TargetsNames)
    println("Initial guess for calibration: ", (exp.(x0) .* logs + x0 .* (.!logs)) .* (10.0 .^ -expon))
    println("Upper bounds for calibration: ", (exp.(upper) .* logs + upper .* (.!logs)) .* (10.0 .^ -expon))
    println("Lower bounds for calibration: ", (exp.(lower) .* logs + lower .* (.!logs)) .* (10.0 .^ -expon))

    PCalib = deepcopy(P)
    OptCalib = deepcopy(Opt)
    @set! OptCalib.display = 0
    @set! OptCalib.tol = 1.0e-7

    G = grids(PCalib, OptCalib)

    function f(x)
        SS = steady_init(G, PCalib, OptCalib)
        crit, flag, _ = calibobjective(x, logs, expon, weights, estvec, Targets, SS, PCalib, OptCalib)
        return crit, flag
    end

    if Opt.leverage_mode == :none
        dir_name = "../temp/no_leverage/"
    elseif Opt.leverage_mode == :riskfree
        dir_name = "../temp/riskfree_leverage/"
    end

    pre_correction_physical = [
        0.016789032387849456,
        1.1760873494390305,
        0.02332592076760511,
        0.41093421143964093,
        2.1013876827525335,
        0.8448572362615245,
        0.1154258911537929,
        0.3015749663377754,
        0.09960006844512176,
        0.12857261301781306,
    ]
    pre_correction_point = similar(pre_correction_physical)
    for ii in eachindex(pre_correction_physical)
        scaled_value = pre_correction_physical[ii] * 10.0^expon[ii]
        pre_correction_point[ii] = logs[ii] ? log(scaled_value) : scaled_value
    end
    initial_points = vcat(reshape(x0, 1, :),
                          reshape(pre_correction_point, 1, :))

    res = MyGlobalOptim.myglobaloptim(x -> f(x), lower, upper, Opt.Ncalib, Opt.Nstarcalib;
        doreset=Opt.resetcalib, doresetinuse=Opt.resetinuse,
        dir_name=dir_name, maxtimeLocal=Opt.maxtimeLocalCalib,
        initial_points=initial_points)

    # Check whether x is within bounds
    for ii in eachindex(res)
        if abs(res[ii] - lower[ii]) < 1e-8
            @warn "Parameter $(ParamNames[ii]) is at lower bound: $(res[ii])"
        elseif abs(res[ii] - upper[ii]) < 1e-8
            @warn "Parameter $(ParamNames[ii]) is at upper bound: $(res[ii])"
        end
    end

    xx = Array{Float64}(undef, length(res))
    xx .= res
    for ii in eachindex(res)
        if logs[ii]
            xx[ii] = exp.(res[ii]) ./ (10.0 .^ expon[ii])
        else
            xx[ii] = res[ii] ./ (10.0 .^ expon[ii])
        end
    end

    x = [PCalib.ρ, PCalib.ρE / PCalib.ρ, PCalib.kappa, PCalib.χ1, PCalib.χ2, PCalib.entre_prod, PCalib.log_y_spread, PCalib.sigma_e, PCalib.χe1, PCalib.χe2, PCalib.alphaaE, PCalib.delta_e + (1 - PCalib.eta / (PCalib.eta + 1)) * PCalib.lambda_e * Opt.entreriskjump, PCalib.lambda_e, PCalib.eta, PCalib.avgtransrate, PCalib.lowstateshare]
    x[estvec] = xx
    #           ρ           ρE        κ             χ1        χ2  entre_prod   log_y_spread    σe          χe1    χe2 αE

    @set! PCalib.ρ = x[1]
    @set! PCalib.ρE = x[2] * x[1]
    @set! PCalib.kappa = x[3]
    @set! PCalib.χ1 = x[4]
    @set! PCalib.χ2 = x[5]
    @set! PCalib.entre_prod = x[6]
    @set! PCalib.log_y_spread = x[7]
    @set! PCalib.sigma_e = x[8]
    @set! PCalib.χe1 = x[9]
    @set! PCalib.χe2 = x[10]
    @set! PCalib.alphaaE = x[11]
    @set! PCalib.lambda_e = x[13]
    @set! PCalib.eta = x[14]
    @set! PCalib.avgtransrate = x[15]
    @set! PCalib.lowstateshare = x[16]
    transition_share_term = PCalib.lowstateshare^2 + (1 - PCalib.lowstateshare)^2
    @set! PCalib.lambda_y_lh = PCalib.avgtransrate * transition_share_term / PCalib.lowstateshare
    @set! PCalib.lambda_y_hl = PCalib.avgtransrate * transition_share_term / (1 - PCalib.lowstateshare)
    @set! PCalib.delta_e = x[12] - (1-PCalib.eta/(PCalib.eta+1))*PCalib.lambda_e*Opt.entreriskjump

    G = grids(PCalib, OptCalib)
    SS, SSall, Pol, PolE, A, Aentre = findsteady(SS, PCalib, G, Opt)
    pHtm_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    wHtm_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    negLiq_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    pHtme_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    wHtme_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    negLiqe_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    simEinc = any(weights[14:18] .> 0.0)
    Stats = calibstats!(SSall, Aentre, PolE, G, PCalib, Opt, pHtm_MAT, wHtm_MAT, negLiq_MAT, pHtme_MAT, wHtme_MAT, negLiqe_MAT, simEinc)

    @save Opt.outpath * "Targets.jld2" Targets

    @save Opt.outpath * "CalibStats.jld2" Stats

    if P.share_entre > 0
        @save Opt.outpath * "ParametersCalib.jld2" PCalib
    else
        @save Opt.outpath * "ParametersNoEntreCalib.jld2" PCalib
    end

    return PCalib, Targets, Stats

end

" `calibobjective(param,logs,expon,W,estvec,Targets, P,Opt)`

`Targets` is a structure containing the calibration targets. `param` is a vector of parameters. `logs` is a vector of the same size as param indicating that the corresponding value in `param` is the log of the actual parameter (this is to be able to impose positive values.)
`W` is a weighting vector for the calibration targets. `estvec` is a vector of booleans indicating which parameters are to be estimated. `SS` is the steady state structure. `P` is the parameter structure. `Opt` is the options structure. The function returns the sum of the squared deviations of the calibration targets from the actual values.
"
function calibobjective(param, logs, expon, weights, estvec, Targets, SSin, P, Opt)
    xx = Array{Float64}(undef, length(param))
    xx .= param
    for ii in eachindex(param)
        if logs[ii]
            xx[ii] = exp.(param[ii]) ./ (10.0 .^ expon[ii])
        else
            xx[ii] = param[ii] ./ (10.0 .^ expon[ii])
        end
    end

    x = [P.ρ, P.ρE / P.ρ, P.kappa, P.χ1, P.χ2, P.entre_prod, P.log_y_spread, P.sigma_e, P.χe1, P.χe2, P.alphaaE, P.delta_e + (1 - P.eta / (P.eta + 1)) * P.lambda_e * Opt.entreriskjump, P.lambda_e, P.eta, P.avgtransrate, P.lowstateshare]
    x[estvec] = xx
    #           ρ           ρE        κ             χ1        χ2  entre_prod   log_y_spread    σe          χe1    χe2 αE

    @set! P.ρ = x[1]
    @set! P.ρE = x[2] * x[1]
    @set! P.kappa = x[3]
    @set! P.χ1 = x[4]
    @set! P.χ2 = x[5]
    @set! P.entre_prod = x[6]
    @set! P.log_y_spread = x[7]
    @set! P.sigma_e = x[8]
    @set! P.χe1 = x[9]
    @set! P.χe2 = x[10]
    @set! P.alphaaE = x[11]
    @set! P.lambda_e = x[13]
    @set! P.eta = x[14]
    @set! P.avgtransrate = x[15]
    @set! P.lowstateshare = x[16]
    transition_share_term = P.lowstateshare^2 + (1 - P.lowstateshare)^2
    @set! P.lambda_y_lh = P.avgtransrate * transition_share_term / P.lowstateshare
    @set! P.lambda_y_hl = P.avgtransrate * transition_share_term / (1 - P.lowstateshare)
    @set! P.delta_e = x[12] - (1-P.eta/(P.eta+1))*P.lambda_e*Opt.entreriskjump

    if P.delta_e < 0
        return Inf, -1.0, zeros(18)
    end

    G = grids(P, Opt) # Rebuild grids for the current calibration parameters.

    pHtm_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    wHtm_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    negLiq_MAT = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    pHtme_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    wHtme_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    negLiqe_MAT = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)

    calibobj = zeros(18)

    try
        if Opt.calibonlyE
            SSall, Aentre, PolE = calibonlyE!(SSin, P, G, Opt)
        else
            _, SSall, _, PolE, _, Aentre = findsteady(SSin, P, G, Opt)
        end

        simEinc = any(weights[14:18] .> 0.0) # whether simulation of entrepreneurial income process is needed
        Statistics = calibstats!(SSall, Aentre, PolE, G, P, Opt, pHtm_MAT, wHtm_MAT, negLiq_MAT, pHtme_MAT, wHtme_MAT, negLiqe_MAT, simEinc)

        calibobj[1] = (Statistics.AggrIlliqSavAnnGDP ./ Targets.AggrIlliqSavAnnGDP .- 1.0)
        calibobj[2] = (Statistics.AggrLiqSavAnnGDP ./ Targets.AggrLiqSavAnnGDP .- 1.0)
        calibobj[3] = (Statistics.pHtmShare ./ Targets.pHtmShare .- 1.0)
        calibobj[4] = (Statistics.wHtmShare ./ Targets.wHtmShare .- 1.0)
        calibobj[5] = Statistics.negLiqShare ./ Targets.negLiqShare .- 1.0
        calibobj[6] = (Statistics.LaborShareE ./ Targets.LaborShareE .- 1.0)
        calibobj[7] = (Statistics.ShareLiqE ./ Targets.ShareLiqE .- 1.0)
        calibobj[8] = (Statistics.ShareIlliqE ./ Targets.ShareIlliqE .- 1.0)
        calibobj[9] = (Statistics.HtmShareE ./ Targets.HtmShareE .- 1.0)
        calibobj[10] = (Statistics.Avg_NW_ratio ./ Targets.Avg_NW_ratio .- 1.0)
        calibobj[11] = (Statistics.AggrCapAnnGDP ./ Targets.AggrCapAnnGDP .- 1.0)
        calibobj[12] = ((Statistics.wHtmShare .+ Statistics.pHtmShare) ./ (Targets.pHtmShare .+ Targets.wHtmShare) .- 1.0)
        calibobj[13] = Statistics.Emp_ratio ./ Targets.Emp_ratio .- 1.0
        calibobj[14] = Statistics.EInc_immobility_ratio ./ Targets.EInc_immobility_ratio .- 1.0
        calibobj[15] = Statistics.EInc_rank_autocorrelation ./ Targets.EInc_rank_autocorrelation .- 1.0
        calibobj[16] = Statistics.EInc_large_loss_probability ./ Targets.EInc_large_loss_probability .- 1.0
        calibobj[17] = Statistics.EInc_large_gain_probability ./ Targets.EInc_large_gain_probability .- 1.0
        calibobj[18] = Statistics.EInc_avg_deciles_lost ./ Targets.EInc_avg_deciles_lost .- 1.0

        crit = calibobj[weights]' * calibobj[weights]

        println("Parameters: ρ = ", P.ρ, ", ρE = ", P.ρE, ", κ = ", P.kappa, ", χ1 = ", P.χ1, ", χ2 = ", P.χ2, ", entre_prod = ", P.entre_prod, ", log_y_spread = ", P.log_y_spread, ", σe = ", P.sigma_e, ", χe1 = ", P.χe1, ", χe2 = ", P.χe2, ", αE = ", P.alphaaE, ", δe = ", P.delta_e, ", λe = ", P.lambda_e, ", η = ", P.eta, ", avgtransrate = ", P.avgtransrate, ", lowstateshare = ", P.lowstateshare, ", λy_lh = ", P.lambda_y_lh, ", λy_hl = ", P.lambda_y_hl)
        printcalib = ["Model" Statistics.AggrIlliqSavAnnGDP Statistics.AggrLiqSavAnnGDP Statistics.pHtmShare Statistics.wHtmShare Statistics.negLiqShare Statistics.LaborShareE Statistics.ShareLiqE Statistics.ShareIlliqE Statistics.HtmShareE Statistics.Avg_NW_ratio Statistics.AggrCapAnnGDP (Statistics.wHtmShare.+Statistics.pHtmShare) Statistics.Emp_ratio Statistics.EInc_immobility_ratio Statistics.EInc_rank_autocorrelation Statistics.EInc_large_loss_probability Statistics.EInc_large_gain_probability Statistics.EInc_avg_deciles_lost;
            "Data" Targets.AggrIlliqSavAnnGDP Targets.AggrLiqSavAnnGDP Targets.pHtmShare Targets.wHtmShare Targets.negLiqShare Targets.LaborShareE Targets.ShareLiqE Targets.ShareIlliqE Targets.HtmShareE Targets.Avg_NW_ratio Targets.AggrCapAnnGDP (Targets.wHtmShare.+Targets.pHtmShare) Targets.Emp_ratio Targets.EInc_immobility_ratio Targets.EInc_rank_autocorrelation Targets.EInc_large_loss_probability Targets.EInc_large_gain_probability Targets.EInc_avg_deciles_lost;
            "Deviation" calibobj[1] calibobj[2] calibobj[3] calibobj[4] calibobj[5] calibobj[6] calibobj[7] calibobj[8] calibobj[9] calibobj[10] calibobj[11] calibobj[12] calibobj[13] calibobj[14] calibobj[15] calibobj[16] calibobj[17] calibobj[18]]
        column_labels = ["", "A/Y", "B/Y", "Phtm", "Whtm", "b < 0", "Lab at e", "Liq e", "Illiq e", "Htm e", "NWratio", "K/Y", "Htmshare", "Emp ratio", "Immobility", "Rank autocorr.", "Pr(loss >= 3 deciles)", "Pr(gain >= 3 deciles)", "Avg. deciles lost | loss >= 3"]
        pretty_table(printcalib, column_labels=column_labels, formatters = [fmt__printf("%5.3f")])

        println("Sum of squared deviations: ", crit)

        return crit, 0.0, calibobj
    catch

        println("Parameters: ρ = ", P.ρ, ", ρE = ", P.ρE, ", κ = ", P.kappa, ", χ1 = ", P.χ1, ", χ2 = ", P.χ2, ", entre_prod = ", P.entre_prod, ", log_y_spread = ", P.log_y_spread, ", σe = ", P.sigma_e, ", χe1 = ", P.χe1, ", χe2 = ", P.χe2, ", αE = ", P.alphaaE, ", δe = ", P.delta_e, ", λe = ", P.lambda_e, ", η = ", P.eta, ", avgtransrate = ", P.avgtransrate, ", lowstateshare = ", P.lowstateshare, ", λy_lh = ", P.lambda_y_lh, ", λy_hl = ", P.lambda_y_hl)
        @warn "Could not find steady state, setting objective to Inf"
        Statistics = false
        crit = Inf
        return crit, -1.0, calibobj
    end
end

# Increment when changes to the bequest equations make earlier saved results
# economically incompatible even if parameters and grids are unchanged.
const BEQUEST_RESULT_VERSION = 1

# Re-solving an unchanged steady state can produce small numerical differences,
# especially when regimes are solved concurrently. Parameters, grids, model
# options are still required to match exactly through the validation signature
# below. Solver changes retain the path but require a new convergence check.
# These tolerances apply only to the freshly
# recomputed terminal steady state.
const BEQUEST_TERMINAL_SCALAR_RTOL = 1e-4
const BEQUEST_TERMINAL_SCALAR_ATOL = 1e-8
const BEQUEST_TERMINAL_DISTRIBUTION_L1_TOL = 1e-3

function bequest_validation_signature(P, G, Opt, beqsize)
    return (
        version=BEQUEST_RESULT_VERSION,
        beqsize=Float64(beqsize),
        parameters=P,
        model_options=(
            entreriskjump=Opt.entreriskjump,
            entreriskdiffusion=Opt.entreriskdiffusion,
            leverage_mode=Opt.leverage_mode,
            stickywages=Opt.stickywages,
            ltau=Opt.ltau,
            a_min2=Opt.a_min2,
        ),
        solver_options=(
            MIT_tol=Opt.MIT_tol,
            bondadjust_MIT=Opt.bondadjust_MIT,
            Entreadjust_MIT=Opt.Entreadjust_MIT,
            stepstickytransK=Opt.stepstickytransK,
            stepstickytransB=Opt.stepstickytransB,
            bicgstab_tol=Opt.bicgstab_tol,
        ),
        grid=(
            N_b=G.N_b,
            N_be=G.N_be,
            N_a=G.N_a,
            N_e=G.N_e,
            N_j=G.N_j,
            N_k=G.N_k,
            N_y=G.N_y,
            N_t=G.N_t,
            b_grid=copy(G.b_grid),
            b_grid_e=copy(G.b_grid_e),
            a_grid=copy(G.a_grid),
            e_grid=copy(G.e_grid),
            y_grid=copy(G.y_grid),
            dt=copy(G.dt),
            mp_shock_t=copy(G.mp_shock_t),
        ),
    )
end

function normalized_l1_distance(x, y)
    size(x) == size(y) || return Inf

    distance = 0.0
    scale = 0.0
    @inbounds for i in eachindex(x, y)
        xi = x[i]
        yi = y[i]
        (isfinite(xi) && isfinite(yi)) || return Inf
        distance += abs(xi - yi)
        scale += max(abs(xi), abs(yi))
    end

    return iszero(scale) ? distance : distance / scale
end

function compatible_bequest_terminal(saved, current)
    scalar_fields = (:rb, :ra, :Output_total, :AggrLiqSav, :AggrIlliqSav,
                     :AggrFirmSav, :BondDem, :zeta_corr)
    scalars_match = all(isapprox(getproperty(saved, field),
                                 getproperty(current, field);
                                 rtol=BEQUEST_TERMINAL_SCALAR_RTOL,
                                 atol=BEQUEST_TERMINAL_SCALAR_ATOL)
                        for field in scalar_fields)
    return scalars_match &&
           normalized_l1_distance(saved.g, current.g) <=
               BEQUEST_TERMINAL_DISTRIBUTION_L1_TOL &&
           normalized_l1_distance(saved.ge, current.ge) <=
               BEQUEST_TERMINAL_DISTRIBUTION_L1_TOL
end

function compatible_bequest_transition(MIT, TransPolW, TransPolE, SSallsaved,
    Psaved, transition_signature, policy_signature, current_signature,
    SSallcurrent, Pcurrent, G)
    return compatible_checkpoint_signature(transition_signature, current_signature) &&
           compatible_checkpoint_signature(policy_signature, current_signature) &&
           isequal(Psaved, Pcurrent) &&
           compatible_bequest_terminal(SSallsaved, SSallcurrent) &&
           length(MIT.rb_t) == G.N_t &&
           length(MIT.mp_shock_t) == G.N_t &&
           length(MIT.AggrLiqSav_t) == G.N_t + 1 &&
           isapprox(MIT.mp_shock_t, G.mp_shock_t; rtol=0.0, atol=0.0) &&
           length(TransPolW) == G.N_t &&
           length(TransPolE) == G.N_t &&
           all(length(policy.c) == G.N for policy in TransPolW) &&
           all(length(policy.ce) == G.Ne for policy in TransPolE)
end

function bequestregimes!(SSallbaseline, Pbeq, Gbeq, Optbeq, beqsize, name,
    AggrLiqSav; n_generators::Integer=0)
    @set! Optbeq.display = 0
    println("Bequest regime: " * string(Pbeq.frac_inh_wealth * 100.0) * "% of entrepreneurs get " * string(beqsize) * "k USD")

    validation_signature =
        bequest_validation_signature(Pbeq, Gbeq, Optbeq, beqsize)
    steady_path = Optbeq.outpath * "SS" * name * ".jld2"
    transition_path = Optbeq.outpath * "MIT_" * name * ".jld2"
    policy_path = Optbeq.outpath * "MIT_policies_" * name * ".jld2"

    loaded_regime_ss = false
    SSallbeqin = deepcopy(SSallbaseline)
    try
        SSallcandidate, saved_signature =
            load(steady_path, "SSallbeq", "validation_signature")
        if compatible_checkpoint_signature(saved_signature, validation_signature)
            SSallbeqin = SSallcandidate
            loaded_regime_ss = true
            println("Loaded validated steady state for bequest regime " * name)
        else
            println("Saved steady state is incompatible with bequest regime " *
                    name * "; starting from baseline steady state")
        end
    catch err
        println("No validated saved steady state found for bequest regime " *
                name * "; starting from baseline steady state (" *
                sprint(showerror, err) * ")")
    end

    SSbeq, SSallbeq, Polbeq, PolEbeq, Abeq, Aentrebeq =
        findsteady_fixedB(SSallbeqin, AggrLiqSav, Pbeq, Gbeq, Optbeq)
    @save steady_path SSbeq SSallbeq Polbeq PolEbeq Abeq Aentrebeq validation_signature

    # Always initialize from the steady state found in this run. In particular,
    # this prevents a baseline-steady-state MIT object from being used when the
    # bequest steady state did not exist on disk.
    MITbeq, TransPolWbeq, TransPolEbeq =
        MIT_init(SSallbeq, Gbeq.mp_shock_t, Polbeq, PolEbeq, Gbeq, Pbeq)

    # A saved transition is only a valid continuation when its corresponding
    # regime steady state was available at the start of this run.
    if loaded_regime_ss
        try
            MITcandidate, SSallsaved, Psaved, transition_signature =
                load(transition_path, "MIT", "SSallend", "P",
                     "validation_signature")
            TransPolWcandidate, TransPolEcandidate, policy_signature =
                load(policy_path, "TransPolW", "TransPolE",
                     "validation_signature")
            if compatible_bequest_transition(
                MITcandidate, TransPolWcandidate, TransPolEcandidate,
                SSallsaved, Psaved, transition_signature, policy_signature,
                validation_signature, SSallbeq, Pbeq, Gbeq)
                MITbeq = MITcandidate
                TransPolWbeq = TransPolWcandidate
                TransPolEbeq = TransPolEcandidate
                println("Loaded validated MIT shock sequence for bequest regime " * name *
                        ", iteration = " * string(MITbeq.it[1]) *
                        ", converged: " * string(MITbeq.converged[1]))
                if !checkpoint_solver_settings_match(transition_signature, validation_signature) ||
                   !checkpoint_solver_settings_match(policy_signature, validation_signature)
                    println("Solver settings changed for bequest regime " * name *
                            "; retaining the saved path and policies and rechecking convergence")
                    MITbeq.it[1] = 0
                    MITbeq.converged[1] = false
                    MITbeq.crit .= Inf
                elseif !MITbeq.converged[1]
                    println("Did not converge for bequest regime " * name * ", resetting iteration counter to 0.")
                    MITbeq.it[1] = 0
                elseif MITbeq.it[1] >= Optbeq.MaxIterMIT
                    println("Validated transition is already at MaxIterMIT; " *
                            "resetting its counter to rebuild generator paths.")
                    MITbeq.it[1] = 0
                    MITbeq.converged[1] = false
                    MITbeq.crit .= Inf
                end
            else
                println("Saved MIT shock sequence or policies are incompatible " *
                        "with bequest regime " * name * "; starting from scratch")
            end
        catch err
            println("No validated MIT shock sequence found for bequest regime " *
                    name * "; starting from scratch (" *
                    sprint(showerror, err) * ")")
        end
    else
        println("Starting MIT shock sequence for bequest regime " * name * " from the newly computed steady state")
    end

    doplot = true
    dosave = true
    entre_exp = 0
    worker_exp = 0
    MaxAbsCapDev = zeros(Optbeq.MaxIterMIT)
    MaxAbsBondDev = zeros(Optbeq.MaxIterMIT)
    @set! Optbeq.display = 1
    TempE = EntreTemp(Gbeq, Optbeq)
    TempW = Temp(Gbeq, Optbeq)
    MITbeq.G_t .= SSallbeq.G

    MITgeneratorsbeq =
        solveMIT!(MITbeq, MaxAbsCapDev, MaxAbsBondDev, TransPolWbeq, TransPolEbeq,
                  SSallbeq.g, SSallbeq.ge, SSallbeq, Gbeq.mp_shock_t, TempW, TempE,
                  false, Gbeq, Pbeq, Optbeq, name, doplot, dosave, entre_exp,
                  worker_exp; n_generators=n_generators,
                  validation_signature=validation_signature)

    MITbeq.converged[1] ||
        error("Bequest transition $name did not converge after " *
              "$(MITbeq.it[1]) iterations; statistics were not computed")

    return (; SSallbeq, Polbeq, PolEbeq, Abeq, Aentrebeq, MITbeq,
            TransPolWbeq, TransPolEbeq, MITgeneratorsbeq)
end

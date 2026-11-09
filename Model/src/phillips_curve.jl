const PHILLIPS_CURVE_RESULT_VERSION = 1

function phillips_curve_validation_signature(P, G, Opt, mp_shock_t, experiment)
    return (
        version=PHILLIPS_CURVE_RESULT_VERSION,
        experiment=String(experiment),
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
            MaxIterMIT=Opt.MaxIterMIT,
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
        ),
        mp_shock_t=copy(mp_shock_t),
    )
end

function compatible_phillips_curve_terminal(saved, current)
    scalar_fields = (:rb, :ra, :Output_total, :AggrLiqSav, :AggrIlliqSav,
                     :AggrFirmSav, :BondDem, :zeta_corr)
    scalars_match = all(isapprox(getproperty(saved, field),
                                 getproperty(current, field);
                                 rtol=1e-8, atol=1e-10)
                        for field in scalar_fields)
    return scalars_match &&
           size(saved.g) == size(current.g) &&
           size(saved.ge) == size(current.ge) &&
           isapprox(saved.g, current.g; rtol=1e-7, atol=1e-8) &&
           isapprox(saved.ge, current.ge; rtol=1e-7, atol=1e-8)
end

function compatible_phillips_curve_transition(MIT, TransPolW, TransPolE,
    SSallsaved, Psaved, transition_signature, policy_signature,
    current_signature, SSallcurrent, Pcurrent, G, mp_shock_t)
    return compatible_checkpoint_signature(transition_signature, current_signature) &&
           compatible_checkpoint_signature(policy_signature, current_signature) &&
           isequal(Psaved, Pcurrent) &&
           compatible_phillips_curve_terminal(SSallsaved, SSallcurrent) &&
           length(MIT.rb_t) == G.N_t &&
           length(MIT.mp_shock_t) == G.N_t &&
           length(MIT.AggrLiqSav_t) == G.N_t + 1 &&
           isapprox(MIT.mp_shock_t, mp_shock_t; rtol=0.0, atol=0.0) &&
           length(TransPolW) == G.N_t &&
           length(TransPolE) == G.N_t &&
           all(length(policy.c) == G.N for policy in TransPolW) &&
           all(length(policy.ce) == G.Ne for policy in TransPolE)
end

function phillips_curve_point_options(Opt, mp_shock_t, experiment)
    OptPoint = Opt
    if mp_shock_t[1] < -0.0025
        @set! OptPoint.stepstickytransK = 0.03
    end

    # shock() converts annual rates to quarterly rates.
    if Opt.leverage_mode == :none && experiment == "bequest_regime_3" &&
       isapprox(mp_shock_t[1], -0.005; rtol=0.0, atol=1e-14)
        @set! OptPoint.stepstickytransB = 0.0005
    end
    return OptPoint
end

# Reuse a compatible checkpoint for this exact shock; otherwise use steady-state guesses.
function solve_phillips_curve_point(SSall, Pol, PolE, mp_shock_t, G, P, Opt,
    name, experiment; doplot::Bool=false, dosave::Bool=true)
    OptPoint = phillips_curve_point_options(Opt, mp_shock_t, experiment)
    if OptPoint.stepstickytransK != Opt.stepstickytransK
        println("Using stepstickytransK = 0.03 for Phillips-curve transition " * name)
    end
    if OptPoint.stepstickytransB != Opt.stepstickytransB
        println("Using stepstickytransB = 0.0005 for Phillips-curve transition " * name)
    end

    validation_signature =
        phillips_curve_validation_signature(P, G, OptPoint, mp_shock_t, experiment)
    transition_path = OptPoint.outpath * "MIT_" * name * ".jld2"
    policy_path = OptPoint.outpath * "MIT_policies_" * name * ".jld2"

    MIT, TransPolW, TransPolE =
        MIT_init(SSall, mp_shock_t, Pol, PolE, G, P)
    MaxAbsCapDev = zeros(OptPoint.MaxIterMIT)
    MaxAbsBondDev = zeros(OptPoint.MaxIterMIT)

    try
        MITcandidate, SSallsaved, Psaved, MaxAbsCapDevcandidate,
        MaxAbsBondDevcandidate, transition_signature =
            load(transition_path, "MIT", "SSallend", "P", "MaxAbsCapDev",
                 "MaxAbsBondDev", "validation_signature")
        TransPolWcandidate, TransPolEcandidate, policy_signature =
            load(policy_path, "TransPolW", "TransPolE", "validation_signature")

        if compatible_phillips_curve_transition(
            MITcandidate, TransPolWcandidate, TransPolEcandidate,
            SSallsaved, Psaved, transition_signature, policy_signature,
            validation_signature, SSall, P, G, mp_shock_t)
            MIT = MITcandidate
            TransPolW = TransPolWcandidate
            TransPolE = TransPolEcandidate
            if length(MaxAbsCapDevcandidate) == OptPoint.MaxIterMIT &&
               length(MaxAbsBondDevcandidate) == OptPoint.MaxIterMIT
                MaxAbsCapDev = MaxAbsCapDevcandidate
                MaxAbsBondDev = MaxAbsBondDevcandidate
            end
            println("Loaded validated Phillips-curve transition " * name *
                    ", iteration = " * string(MIT.it[1]) *
                    ", converged: " * string(MIT.converged[1]))
            if !checkpoint_solver_settings_match(transition_signature, validation_signature) ||
               !checkpoint_solver_settings_match(policy_signature, validation_signature)
                println("Solver settings changed for " * name *
                        "; retaining the saved path and policies and rechecking convergence")
                MIT.it[1] = 0
                MIT.converged[1] = false
                MIT.crit .= Inf
                MaxAbsCapDev .= 0
                MaxAbsBondDev .= 0
            end
        else
            println("Saved Phillips-curve transition " * name *
                    " is incompatible; starting from its steady state")
        end
    catch err
        println("No validated saved Phillips-curve transition found for " *
                name * "; starting from its steady state (" *
                sprint(showerror, err) * ")")
    end

    if MIT.converged[1]
        return (; MIT, MaxAbsCapDev, MaxAbsBondDev)
    end

    if MIT.it[1] > 0
        println("Phillips-curve transition " * name *
                " is a validated but unconverged warm start; resetting its iteration counter")
        MIT.it[1] = 0
        MIT.converged[1] = false
        MIT.crit .= Inf
    end

    entre_exp = 0
    worker_exp = 0
    MIT, TransPolW, TransPolE, MaxAbsCapDev, MaxAbsBondDev =
        solveMIT(MIT, TransPolW, TransPolE, SSall.g, SSall.ge, SSall,
                 mp_shock_t, false, G, P, OptPoint, name, doplot, dosave,
                 entre_exp, worker_exp; copy_inputs=false,
                 validation_signature=validation_signature)

    return (; MIT, MaxAbsCapDev, MaxAbsBondDev)
end

function solve_phillips_curve_grid(SSall, Pol, PolE, G, P, Opt;
    shocks=collect(-0.02:0.005:0.02), name_prefix="", experiment,
    output_path, doplot::Bool=false, dosave::Bool=true)
    mp_shock_tVec = [shock(sval, exp(-0.5), G.N_t, G.dt) for sval in shocks]
    namesVec = [name_prefix * "shocksize" * string(path[1])
                for path in mp_shock_tVec]

    solve_point = (mp_shock_t, name) ->
        solve_phillips_curve_point(SSall, Pol, PolE, mp_shock_t, G, P, Opt,
                                   name, experiment; doplot=doplot,
                                   dosave=dosave)
    results = nprocs() > 1 ?
        pmap(solve_point, mp_shock_tVec, namesVec) :
        map(solve_point, mp_shock_tVec, namesVec)

    MITVec = [result.MIT for result in results]
    MaxAbsCapDevVec = [result.MaxAbsCapDev for result in results]
    MaxAbsBondDevVec = [result.MaxAbsBondDev for result in results]

    failed = findall(mit -> !mit.converged[1], MITVec)
    isempty(failed) ||
        error("Phillips-curve transitions did not converge at shock indices " *
              join(failed, ", ") * "; aggregate results were not saved")

    validation_signatures =
        [phillips_curve_validation_signature(P, G,
            phillips_curve_point_options(Opt, mp_shock_t, experiment), mp_shock_t, experiment)
         for mp_shock_t in mp_shock_tVec]
    @save output_path MITVec MaxAbsCapDevVec MaxAbsBondDevVec shocks mp_shock_tVec namesVec validation_signatures
    return MITVec, MaxAbsCapDevVec, MaxAbsBondDevVec
end

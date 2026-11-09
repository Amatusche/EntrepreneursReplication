function prepare_policy_inputs(TransPolWInitGuess, TransPolEInitGuess, copy_inputs::Bool)
    if copy_inputs
        return deepcopy(TransPolWInitGuess), deepcopy(TransPolEInitGuess)
    end
    return TransPolWInitGuess, TransPolEInitGuess
end

"`solveMIT(MITInitGuess,TransPolWInitGuess,TransPolEInitGuess,gstart,gestart,SSallend,mp_shock_t,A0fix,G,P,Opt,ModelName,doplot,dosave,entre_exp,worker_exp; copy_inputs=true)`
Solves MIT shock from initial distributions `gstart` and `gestart` with final steady state `SSallend`.
Initializes `MIT` and convergence criteria (`MaxAbsCapDev`,`MaxAbsBondDev`) .
`A0fix` = true to keep initial distributions unaffected by stock price q.
Set `copy_inputs=false` to update the supplied policy paths in place and avoid
duplicating them in memory. The caller must not rely on the original policy
values after such a call. `MITInitGuess` is always copied."
function solveMIT(MITInitGuess,TransPolWInitGuess,TransPolEInitGuess,gstart,gestart,SSallend,mp_shock_t,A0fix,G,P,Opt,ModelName,doplot,dosave,entre_exp,worker_exp,gstartNA=false,shareNA=0.0; n_generators::Integer=0, copy_inputs::Bool=true, validation_signature=nothing)

    MIT = deepcopy(MITInitGuess)
    TransPolW, TransPolE =
        prepare_policy_inputs(TransPolWInitGuess, TransPolEInitGuess, copy_inputs)
    TempE = EntreTemp(G, Opt)
    TempW = Temp(G, Opt)
    MaxAbsCapDev = zeros(Opt.MaxIterMIT); # this initializes the convergence criterion for savings
    MaxAbsBondDev = zeros(Opt.MaxIterMIT); # this initializes the convergence criterion for labor supply
    if worker_exp == 4
        TransPolWNA = deepcopy(TransPolWInitGuess)
        generators = solveMIT!(MIT,MaxAbsCapDev,MaxAbsBondDev,TransPolW,TransPolE,gstart,gestart,SSallend,mp_shock_t,TempW,TempE,A0fix,G,P,Opt,ModelName,doplot,dosave,entre_exp,worker_exp,TransPolWNA,gstartNA,shareNA; n_generators=n_generators, validation_signature=validation_signature)
    else
        generators = solveMIT!(MIT,MaxAbsCapDev,MaxAbsBondDev,TransPolW,TransPolE,gstart,gestart,SSallend,mp_shock_t,TempW,TempE,A0fix,G,P,Opt,ModelName,doplot,dosave,entre_exp,worker_exp; n_generators=n_generators, validation_signature=validation_signature)
    end

    if n_generators == 0
        return MIT, TransPolW, TransPolE, MaxAbsCapDev, MaxAbsBondDev
    end
    return MIT, TransPolW, TransPolE, MaxAbsCapDev, MaxAbsBondDev, generators
end

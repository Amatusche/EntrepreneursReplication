function options()

    entreriskjump = false
    leverage_mode = :none # :riskfree, :risky, :none; risk-free debt does not amplify firm-equity risk
    @assert leverage_mode in (:none, :riskfree, :risky)
    entreriskdiffusion = true

    Ncalib = 100000
    Nstarcalib = 20
    maxitersLocalCalib = 150
    reltolLocalCalib = 1e-4
    xtolrelLocalCalib = 1e-3
    resetcalib = false
    resetinuse = false
    resetfailedcalib = true # retry local searches that previously hit a limit; preserves Sobol evaluations
    debacker_seed = 1234 # fixed seed for deterministic entrepreneurial-income simulations

    parallel = false
    stickywages = false
    calibrate = false
    calibonlyE = false # solve only entrepreneur problem for calibration (no GE!)
    RANK = false

    if stickywages
        outpath = "../out/stickywages/"
    else
        outpath = "../out/stickyprices/"
    end
    if RANK
        outpath = outpath*"RANK/"
    end
    if leverage_mode == :none
        outpath = outpath*"no_leverage/"
    elseif leverage_mode == :riskfree
        outpath = outpath*"riskfree_leverage/"
    elseif leverage_mode == :risky
        outpath = outpath*"risky_leverage/"
    end

    sep_inc_states = 0 # set to 1 to update value function separately for different income states (use explicit method in income dimensions and implicit method in asset dimensions)
    if sep_inc_states == 0
        Delta = 1000.0
    elseif sep_inc_states == 1
        Delta = 1.0
    end
    sep_inc_states_e = 0 # set to 1 to update value function separately for different income states (use explicit method in income dimensions and implicit method in asset dimensions)
    if sep_inc_states_e == 0
        Delta_e = 1000.0
    elseif sep_inc_states_e == 1
        Delta_e = 1.0
    end

    do_first_set_only = 0 # only calibrate the five HANK parameters (rho,kappa,chi0,chi1,chi2)
    do_second_set_only = 0 # only calibrate entrepreneurial parameters (entre_prod,nu)

Opt = (
    entreriskjump = entreriskjump,
    entreriskdiffusion = entreriskdiffusion,
    a_min2 = 0.03*2.92/4, # HANK calibration
    sep_inc_states = sep_inc_states, # set to 1 to update value function separately for different income states (use explicit method in income dimensions and implicit method in asset dimensions)
    sep_inc_states_e = sep_inc_states_e,
    Delta_e = Delta_e,
    Delta = Delta,
    Vtol = 1e-10,
    tol = 1e-9,
    # :bicgstabl, :bicgstab, :dqgmres, :conjugategradients, or :direct
    worker_HJB_solver = :bicgstab,
    worker_KFE_solver = :bicgstab,
    entrepreneur_HJB_solver = :bicgstab,
    entrepreneur_KFE_solver = :bicgstab,
    bicgstab_tol = 1e-10,
    MaxIterHJB = 40,
    MaxIterSteady = 50,
    min_c = 1e-6,
    max_d = 1e8,
    v_a_min = 1e-8,
    v_b_min = 1e-8,
    v_e_min = 1e-8,
    caplab_stepsize = 0.05,

    display = 1, # 2 prints HJB iteration details

    use_solver = 1,
    ltau = 18,

    do_first_set_only = do_first_set_only, # only calibrate the five HANK parameters (rho,kappa,chi0,chi1,chi2)
    do_second_set_only = do_first_set_only, # only calibrate entrepreneurial parameters (entre_prod,nu)
    do_all_targets = 1 - do_first_set_only - do_second_set_only, # calibrate all parameters
    do_inh_targets = 0,

    N_DeBackerSim1 = 1_000_000, # Number of households simulated in one matrix
    N_DeBackerSim2= 1, # Number of households matrices simulated, -> Total sample size is N_DeBackerSim*N_DeBackerSimRep
    dt_DeBackerSim = 0.01, # Length of simulation time steps as fraction / multiple of quarters
    debacker_seed = debacker_seed,

    # regarding MIT shock
    EntreSteadyPol = 0,
    EntreSteadyrb = 0,
    bondadjust_MIT = 0.1,
    Entreadjust_MIT = 0.1,
    BiCGSTAB_MIT = 1,
    MIT_tol = 5e-6,
    MaxIterMIT = 500,
    stepstickytransK  = 0.025,
    stepstickytransB  = 0.001,
    parallel = parallel,
    stickywages = stickywages,
    outpath = outpath,
    calibrate = calibrate,
    calibonlyE = calibonlyE, # solve only entrepreneur problem for calibration (no GE!)
    Ncalib = Ncalib,
    Nstarcalib = Nstarcalib,
    maxitersLocalCalib = maxitersLocalCalib,
    reltolLocalCalib = reltolLocalCalib,
    xtolrelLocalCalib = xtolrelLocalCalib,
    resetcalib = resetcalib,
    resetfailedcalib = resetfailedcalib,
    maxtimeLocalCalib = 60.0*60.0*24.0, # maximum time in seconds for local optimization in calibration
    resetinuse = resetinuse,
    leverage_mode = leverage_mode, # set to :riskfree, :risky, or :none to use leverage in the model (firm profits go into liquid account)
)

end

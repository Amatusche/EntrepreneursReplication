using ILUZero, ConjugateGradients

const SteadyKrylovWorkspace = Union{
    Nothing,
    Krylov.BicgstabWorkspace{Float64,Float64,Vector{Float64}},
    Krylov.DqgmresWorkspace{Float64,Float64,Vector{Float64}},
}
const LinearILU = Union{
    Nothing,
    ILUZero.ILU0Precon{Float64,Int64,Float64},
}

# Construct the selected workspace once as part of Temp or EntreTemp.
function steady_krylov_workspace(solver::Symbol, b)
    if solver === :bicgstab
        return Krylov.krylov_workspace(Val(:bicgstab), Krylov.KrylovConstructor(b))
    elseif solver === :dqgmres
        return Krylov.krylov_workspace(
            Val(:dqgmres), Krylov.KrylovConstructor(b); memory=20)
    elseif solver === :bicgstabl || solver === :conjugategradients ||
           solver === :direct
        return nothing
    end
    throw(ArgumentError(
        "Unknown solver $solver; use :bicgstabl, :bicgstab, :dqgmres, " *
        ":conjugategradients, or :direct"))
end

struct SSminstruct
    v::Array{Float64, 4}
    ve::Array{Float64, 4}
    lab0::Array{Float64, 4}
    TransW::Array{Float64, 4}
    TransE::Array{Float64, 4}
    gtemp::Array{Float64, 1}
    gtempE::Array{Float64, 1}
    rb::Float64
    ra::Float64
    pb::Float64
    pe::Float64
    pp::Float64
    w::Float64
    what::Float64
    Output_total::Float64
    AggrOutputE::Float64
    RealProfits::Float64
    AggrCapSup::Float64
    CapLabRatio::Float64
    omega_entre::Float64
    zeta_corr::Float64
end

struct SSstruct
    v::Array{Float64, 4}
    ve::Array{Float64, 4}
    lab0::Array{Float64, 4}
    TransW::Array{Float64, 4}
    TransE::Array{Float64, 4}
    Trans::Float64
    gtemp::Array{Float64, 1}
    gtempE::Array{Float64, 1}
    rb::Float64
    ra::Float64
    pb::Float64
    pe::Float64
    pp::Float64
    w::Float64
    what::Float64
    Output_total::Float64
    AggrOutputE::Float64
    RealProfits::Float64
    AggrCapSup::Float64 # total market-intermediated capital stock supplied through illiquid assets
    KW::Float64 # market-intermediated capital held by workers -> attributed by illiquid-savings shares
    KE::Float64 # market-intermediated capital held by entrepreneurs -> attributed by illiquid-savings shares
    CapLabRatio::Float64
    omega_entre::Float64
    zeta_corr::Float64
    LabDemE::Array{Float64, 4}
    AggrCapDemE::Float64 # market-intermediated capital used in entrepreneur firms, i.e. leveraged (θ-1)e
    AggrCapDemC::Float64 # market-intermediated capital used in the corporate sector
    OutputE::Array{Float64, 4}
    RealGrossProfitsE::Array{Float64, 4}
    TaxableBusinessIncomeE::Array{Float64, 4} # business income before tax, net of expected depreciation
    g::Array{Float64, 4}
    ge ::Array{Float64, 4}
    ge_marg_b::Array{Float64, 1}
    ge_marg_a::Array{Float64, 1}
    ge_marg_e::Array{Float64, 1}
    ge_marg_ey::Array{Float64, 2}
    ge_marg_ba::Array{Float64, 2}
    ge_marg_bae::Array{Float64, 3}
    g_marg_b::Array{Float64, 1}
    g_marg_a::Array{Float64, 1}
    g_marg_ba::Array{Float64, 2}
    AggrLiqSavW::Float64
    AggrIlliqSavW::Float64
    AggrLiqSavE::Float64
    AggrIlliqSavE::Float64
    AggrFirmSav::Float64 # entrepreneurs' own private business capital stock e
    AggrLabDemE::Float64
    AggrLabSup::Float64
    AggrLiqSav::Float64
    AggrIlliqSav::Float64
    AggrLabDemC::Float64
    rk::Float64
    OutputC::Float64
    AdjCosts::Float64
    PortfAdjW::Float64
    PortfAdjE::Float64
    CapAdjE::Float64
    BorrCostsW::Float64
    BorrCostsE::Float64
    BorrCosts::Float64
    PortfAdj::Float64
    BondDem::Float64
    q::Float64
    tax_rev_entre::Float64
    tax_rev_worker::Float64
    subs_payment_union::Float64
    G::Float64
    Inv::Float64 # total investment
    Inv_M::Float64 # total investment in market-intermediated capital
    InvC_M::Float64 # market-intermediated investment used in the corporate sector
    InvE_M::Float64 # market-intermediated investment used in the entrepreneurial sector
    InvE_nonM::Float64 # non-market investment in entrepreneurs' own private business capital
    InvCW::Float64 # worker-owned share of market-intermediated investment
    InvCE::Float64 # entrepreneur-owned share of market-intermediated investment
    InvE::Float64 # total investment used in the entrepreneurial sector, market and non-market
    ConsumpE ::Float64
    ConsumpW::Float64
    C::Float64
    RealProfitsToW::Float64
    RealProfitsToE::Float64
    AvgPerCapitaRealProfitsToW::Float64
    AvgPerCapitaRealProfitsToE::Float64
    UnionProfits::Float64 # Union profits
    Tunion::Float64 # Union transfers
end

struct SteadyStats
    AggrIlliqSavAnnGDP::Float64
    AggrLiqSavAnnGDP::Float64
    pHtmShare::Float64
    wHtmShare::Float64
    negLiqShare::Float64
    LaborShareE::Float64
    PrivateFirmCapShare::Float64
    ShareLiqE::Float64
    ShareIlliqE::Float64
    HtmShareE::Float64
    MedReturn::Array{Float64,2} # taxable business income over business equity
    MedReturnCF::Array{Float64,2} # effective owner cash flow over business equity
    ShareEByNWGroup::Array{Float64,2}
    NWShares::Array{Float64,2}
    DeBackerModelQuint::Array{Float64,2}
    DeBackerModelFull::Array{Float64,2}
    Avg_NW_ratio::Float64
    Emp_ratio::Float64
end

struct Polstruct
    c::Array{Float64, 1}
    lab::Array{Float64, 1}
    d::Array{Float64, 1}
end

struct PolEstruct
    ce::Array{Float64, 1}
    f::Array{Float64, 1}
    de::Array{Float64, 1}
    θ::Array{Float64, 1} # leverage
end

mutable struct EntreTemp

    ce_B::Array{Float64,4}
    ce_F::Array{Float64,4}

    θ_B::Array{Float64,4}
    θ_F::Array{Float64,4}
    θ_0::Array{Float64,4}
    θnew::Array{Float64,4}

    b::Vector{Float64}
    ve_new::Vector{Float64}

    dBB::Array{Float64,4}
    dBF::Array{Float64,4}
    dFB::Array{Float64,4}
    dFF::Array{Float64,4}

    fBB::Array{Float64,4}
    fBF::Array{Float64,4}
    fFB::Array{Float64,4}
    fFF::Array{Float64,4}

    d_B::Array{Float64,4}
    d_F::Array{Float64,4}
    f_B::Array{Float64,4}
    f_F::Array{Float64,4}

    sf_B::Array{Float64,4}
    sf_F::Array{Float64,4}
    sce_B::Array{Float64,4}
    sce_F::Array{Float64,4}
    sde_B::Array{Float64,4}
    sde_F::Array{Float64,4}

    Ic_F::BitArray{4}
    Ic_B::BitArray{4}
    Ic_0::BitArray{4}
    Id_F::BitArray{4}
    Id_B::BitArray{4}
    If_F::BitArray{4}
    If_B::BitArray{4}

    dB::Array{Float64,4}
    dF::Array{Float64,4}
    fB::Array{Float64,4}
    fF::Array{Float64,4}

    ve_b_F::Array{Float64,4}
    ve_b_B::Array{Float64,4}
    ve_a_F::Array{Float64,4}
    ve_a_B::Array{Float64,4}
    ve_e_F::Array{Float64,4}
    ve_e_B::Array{Float64,4}
    ve_e_F_raw::Array{Float64,4}
    ve_e_B_raw::Array{Float64,4}
    ve_ee::Array{Float64,4}

    edriftB::Array{Float64,4}
    edriftF::Array{Float64,4}
    edrift0::Array{Float64,4}

    # Workspaces used by transit_mat_entre!
    trans_Xb::Vector{Float64}
    trans_Y::Vector{Float64}
    trans_Zb::Vector{Float64}
    trans_Xa::Vector{Float64}
    trans_Za::Vector{Float64}
    trans_Xe::Vector{Float64}
    trans_Ze::Vector{Float64}
    trans_B::Array{Float64,4}
    trans_F::Array{Float64,4}

    BCGD::ConjugateGradients.BiCGStabData{Float64}
    HJBILU::LinearILU
    KFEILU::LinearILU
    HJBKrylovWork::SteadyKrylovWorkspace
    KFEKrylovWork::SteadyKrylovWorkspace
    KFEKrylovShadow::Vector{Float64}
end

function EntreTemp(G, Opt)
    Ne = G.Ne
    dims4 = (G.N_be, G.N_a, G.N_e, G.N_y)
    krylov_template = Vector{Float64}(undef, Ne)
    EntreTemp(

        zeros(dims4...),
        zeros(dims4...),

        ones(dims4...),
        ones(dims4...),
        ones(dims4...),
        ones(dims4...),

        zeros(Ne),
        zeros(Ne),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        BitArray{4}(undef, dims4...),
        BitArray{4}(undef, dims4...),
        BitArray{4}(undef, dims4...),
        BitArray{4}(undef, dims4...),
        BitArray{4}(undef, dims4...),
        BitArray{4}(undef, dims4...),
        BitArray{4}(undef, dims4...),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        zeros(dims4...),
        zeros(dims4...),
        zeros(dims4...),

        Vector{Float64}(undef, Ne - 1),
        Vector{Float64}(undef, Ne),
        Vector{Float64}(undef, Ne - 1),
        Vector{Float64}(undef, Ne - G.N_be),
        Vector{Float64}(undef, Ne - G.N_be),
        Vector{Float64}(undef, Ne - G.N_abe),
        Vector{Float64}(undef, Ne - G.N_abe),
        Array{Float64}(undef, dims4...),
        Array{Float64}(undef, dims4...),

        ConjugateGradients.BiCGStabData(Ne,Float64),
        nothing,
        nothing,
        steady_krylov_workspace(Opt.entrepreneur_HJB_solver, krylov_template),
        steady_krylov_workspace(Opt.entrepreneur_KFE_solver, krylov_template),
        Float64[]
    )
end

mutable struct Temp

    c_0::Array{Float64,4}
    c_F::Array{Float64,4}
    c_B::Array{Float64,4}

    lab_F::Array{Float64,4}
    lab_B::Array{Float64,4}
    lab::Array{Float64,4}

    b::Array{Float64,1}
    v_new::Array{Float64,1}

    dBB::Array{Float64,4}
    dBF::Array{Float64,4}
    dFB::Array{Float64,4}
    dFF::Array{Float64,4}

    d_B::Array{Float64,4}
    d_F::Array{Float64,4}

    sc_B::Array{Float64,4}
    sc_F::Array{Float64,4}
    sd_B::Array{Float64,4}
    sd_F::Array{Float64,4}

    Ic_F::BitArray{4}
    Ic_B::BitArray{4}
    Ic_0::BitArray{4}
    Id_F::BitArray{4}
    Id_B::BitArray{4}

    dB::Array{Float64,4}
    dF::Array{Float64,4}

    v_prime_b_F::Array{Float64,4}
    v_prime_b_B::Array{Float64,4}
    v_prime_a_F::Array{Float64,4}
    v_prime_a_B::Array{Float64,4}

    # Workspaces used by transit_mat!
    trans_Xb::Vector{Float64}
    trans_Y::Vector{Float64}
    trans_Zb::Vector{Float64}
    trans_Xa::Vector{Float64}
    trans_Za::Vector{Float64}
    trans_B::Array{Float64,4}
    trans_F::Array{Float64,4}

    BCGD::ConjugateGradients.BiCGStabData{Float64}
    HJBILU::LinearILU
    KFEILU::LinearILU
    HJBKrylovWork::SteadyKrylovWorkspace
    KFEKrylovWork::SteadyKrylovWorkspace
    KFEKrylovShadow::Vector{Float64}

end

function Temp(G, Opt)

    N = G.N
    dims4 = (G.N_b, G.N_a, G.N_j, G.N_k)
    krylov_template = Vector{Float64}(undef, N)

    Temp(

    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),

    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),

    Array{Float64}(undef, N),
    Array{Float64}(undef, N),

    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),

    zeros(dims4...),
    zeros(dims4...),

    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),

    BitArray{4}(undef, dims4...),
    BitArray{4}(undef, dims4...),
    BitArray{4}(undef, dims4...),
    BitArray{4}(undef, dims4...),
    BitArray{4}(undef, dims4...),

    zeros(dims4...),
    zeros(dims4...),

    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),
    zeros(dims4...),

    Vector{Float64}(undef, N - 1),
    Vector{Float64}(undef, N),
    Vector{Float64}(undef, N - 1),
    Vector{Float64}(undef, N - G.N_b),
    Vector{Float64}(undef, N - G.N_b),
    Array{Float64}(undef, dims4...),
    Array{Float64}(undef, dims4...),

    ConjugateGradients.BiCGStabData(N,Float64),
    nothing,
    nothing,
    steady_krylov_workspace(Opt.worker_HJB_solver, krylov_template),
    steady_krylov_workspace(Opt.worker_KFE_solver, krylov_template),
    Float64[]

    )

end

struct MITstruct
    mp_shock_t::Array{Float64, 1}
    rb_t::Array{Float64, 1}
    q_t::Array{Float64, 1}
    mc_t::Array{Float64, 1}
    wagemarkdown_t::Array{Float64, 1}
    rk_t::Array{Float64, 1}
    AggrCapSup_t::Array{Float64, 1} # total market-intermediated capital stock supplied through illiquid assets
    KW_t::Array{Float64, 1} # worker-owned share of market-intermediated capital
    KE_t::Array{Float64, 1} # entrepreneur-owned share of market-intermediated capital
    Inv_t::Array{Float64, 1}
    Inv_M_t::Array{Float64, 1} # total market-intermediated investment
    InvE_t::Array{Float64, 1} # total investment used in the entrepreneurial sector, market and non-market
    InvE_nonM_t::Array{Float64, 1} # non-market investment in entrepreneurs' own private business capital
    InvE_M_t::Array{Float64, 1} # market-intermediated investment used in the entrepreneurial sector
    InvC_M_t::Array{Float64, 1} # market-intermediated investment used in the corporate sector
    InvCW_t::Array{Float64, 1} # worker-owned share of market-intermediated investment
    InvCE_t::Array{Float64, 1} # entrepreneur-owned share of market-intermediated investment
    infl_t::Array{Float64, 1}
    inflw_t::Array{Float64, 1}
    AggrCapDemC_t::Array{Float64, 1} # market-intermediated capital used in the corporate sector
    AggrCapDemE_t::Array{Float64, 1} # market-intermediated capital used in entrepreneur firms, i.e. leveraged (θ-1)e
    AggrLabDemC_t::Array{Float64, 1}
    AggrLabDemE_t::Array{Float64, 1}
    AggrLabSup_t::Array{Float64, 1}
    AggrLiqSav_t::Array{Float64, 1}
    AggrLiqSavW_t::Array{Float64, 1}
    AggrLiqSavE_t::Array{Float64, 1}
    AggrIlliqSav_t::Array{Float64, 1}
    AggrIlliqSavW_t::Array{Float64, 1}
    AggrIlliqSavE_t::Array{Float64, 1}

    ConsumpW_t::Array{Float64, 1}
    ConsumpE_t::Array{Float64, 1}
    PortfAdjW_t::Array{Float64, 1}
    PortfAdjE_t::Array{Float64, 1}
    CapAdjE_t::Array{Float64, 1}
    BorrCosts_t::Array{Float64, 1}
    PortfAdj_t::Array{Float64, 1}

    AggrOutputE_t::Array{Float64, 1}
    AggrFirmSav_t::Array{Float64, 1} # entrepreneurs' own private business capital stock e
    tax_rev_entre_t::Array{Float64, 1}
    subs_payment_union_t::Array{Float64, 1}
    Tunion_t::Array{Float64, 1}

    OutputC_t::Array{Float64, 1}
    Output_total_t::Array{Float64, 1}
    omega_entre::Array{Float64, 1}
    zeta_corr_t::Array{Float64, 1}
    dinfl_t::Array{Float64, 1}
    dOutput_t::Array{Float64, 1}
    dinflw_t::Array{Float64, 1}
    dLab_t::Array{Float64, 1}
    pe_t::Array{Float64, 1}
    w_t::Array{Float64, 1}
    what_t::Array{Float64, 1}
    AdjCosts_t::Array{Float64, 1}
    RealProfits_t::Array{Float64, 1}
    ra_t::Array{Float64, 1}
    Trans_t::Array{Float64, 1}
    G_t::Array{Float64, 1}
    C_t::Array{Float64, 1}
    Excess_Sup_Liq_t::Array{Float64, 1}

    RealProfitsToW_t::Array{Float64, 1}
    RealProfitsToE_t::Array{Float64, 1}
    AvgPerCapitaRealProfitsToW_t::Array{Float64, 1}
    AvgPerCapitaRealProfitsToE_t::Array{Float64, 1}
    UnionProfits_t::Array{Float64, 1}

    it::Array{Int64, 1}
    converged::Array{Bool, 1}
    crit::Array{Float64, 1}
end

struct CalibTargets
    AggrIlliqSavAnnGDP::Float64
    AggrLiqSavAnnGDP::Float64
    pHtmShare::Float64
    wHtmShare::Float64
    negLiqShare::Float64

    LaborShareE::Float64

    ShareLiqE::Float64
    ShareIlliqE::Float64
    HtmShareE::Float64
    Avg_NW_ratio::Float64
    AggrCapAnnGDP::Float64
    Emp_ratio::Float64
    EInc_immobility_ratio::Float64
    EInc_rank_autocorrelation::Float64
    EInc_large_loss_probability::Float64
    EInc_large_gain_probability::Float64
    EInc_avg_deciles_lost::Float64
end

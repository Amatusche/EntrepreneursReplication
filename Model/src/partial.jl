" MIT_partial = partial(MIT, TransPolW, TransPolE, SSall, prices, P, G, Opt)

This function computes partial equilibrium effects of a change of the prices specified in `prices`.

"
function partial(MIT, TransPolW, TransPolE, SSall, prices, price_seq, P, G, Opt; n_generators::Integer=0)
    println("Computing partial equilibrium responses")

    # Initialize

    g = Array{Float64}(undef, G.N_b, G.N_a, G.N_j, G.N_k)
    ge = Array{Float64}(undef, G.N_be, G.N_a, G.N_e, G.N_y)
    AggrOutputE_new = Array{Float64}(undef, G.N_t)  # aggregate output of entrepreneurs
    AggrLabDemE_new = Array{Float64}(undef, G.N_t)  # aggregate labor demand of entrepreneurs
    tax_rev_entre_new = Array{Float64}(undef, G.N_t)
    g_tilde = Array{Float64}(undef, G.N)
    g_tilde .= SSall.g[:]
    g_tilde_e = Array{Float64}(undef, G.Ne)
    g_tilde_e .= SSall.ge[:]
    TempW = Temp(G, Opt)
    TempE = EntreTemp(G, Opt)

    Ainit = deepcopy(G.Aexo)
    # The generators have fixed CSC structure; only their nzval arrays differ by date.
    A_t = sparse_value_path(G.Az, G.N_t)
    Aentre_local_t = sparse_value_path(Ainit, G.N_t)

    MIT_partial = deepcopy(MIT)

    PartialPolW = deepcopy(TransPolW)
    PartialPolE = deepcopy(TransPolE)

    length(prices) == length(price_seq) ||
        throw(ArgumentError("prices and price_seq must have the same length"))

    # Start from steady-state paths and overwrite only the prices included in
    # the partial-equilibrium experiment.
    MIT_partial.rb_t .= SSall.rb
    MIT_partial.w_t .= SSall.w
    MIT_partial.what_t .= SSall.what
    MIT_partial.ra_t .= SSall.ra
    MIT_partial.rk_t .= SSall.rk
    MIT_partial.pe_t .= SSall.pe
    MIT_partial.Trans_t .= SSall.Trans
    MIT_partial.RealProfits_t .= SSall.RealProfits
    MIT_partial.q_t .= SSall.q

    MIT_partial.omega_entre .= SSall.omega_entre
    MIT_partial.AggrOutputE_t .= SSall.AggrOutputE
    MIT_partial.UnionProfits_t .= SSall.UnionProfits
    MIT_partial.zeta_corr_t .= SSall.zeta_corr

    for pp = eachindex(prices)
        if prices[pp] == "rb"
            MIT_partial.rb_t .= price_seq[pp]
        elseif prices[pp] == "w"
            MIT_partial.w_t .= price_seq[pp]
        elseif prices[pp] == "what"
            MIT_partial.what_t .= price_seq[pp]
        elseif prices[pp] == "ra"
            MIT_partial.ra_t .= price_seq[pp]
        elseif prices[pp] == "rk"
            MIT_partial.rk_t .= price_seq[pp]
        elseif prices[pp] == "pe"
            MIT_partial.pe_t .= price_seq[pp]
        elseif prices[pp] == "T_gov"
            MIT_partial.Trans_t .= price_seq[pp]
        elseif prices[pp] == "Profits"
            MIT_partial.RealProfits_t .= price_seq[pp]
        elseif prices[pp] == "q"
            MIT_partial.q_t .= price_seq[pp]
        else
            throw(ArgumentError("unknown partial-equilibrium price: $(prices[pp])"))
        end
    end

    # Initial physical capital is unchanged; only the selected q path can
    # revalue initial illiquid wealth. Do not inherit initial GE asset stocks.
    g .= SSall.g
    ge .= SSall.ge
    A_0_drop = 1.0 + (MIT_partial.q_t[1] - SSall.q) / SSall.AggrIlliqSav
    if A_0_drop != 1.0
        initial_q_change!(g, A_0_drop, G.a_grid, G.N_a, 2)
        if P.share_entre > 0.0
            initial_q_change!(ge, A_0_drop, G.a_grid, G.N_a, 2)
        end
    end
    g_tilde .= vec(g)
    g_tilde_e .= vec(ge)
    MIT_partial.AggrIlliqSavW_t[1] = (1.0 - P.share_entre) * dot(G.a_grid_MAT, g)
    MIT_partial.AggrIlliqSavE_t[1] = P.share_entre * dot(G.a_grid_MAT_e, ge)
    MIT_partial.AggrIlliqSav_t[1] = MIT_partial.AggrIlliqSavW_t[1] + MIT_partial.AggrIlliqSavE_t[1]
    MIT_partial.AggrLiqSavW_t[1] = (1.0 - P.share_entre) * dot(G.b_grid_MAT, g)
    MIT_partial.AggrLiqSavE_t[1] = P.share_entre * dot(G.b_grid_MAT_e, ge)
    MIT_partial.AggrLiqSav_t[1] = MIT_partial.AggrLiqSavW_t[1] + MIT_partial.AggrLiqSavE_t[1]
    MIT_partial.AggrFirmSav_t[1] = P.share_entre * dot(G.e_grid_MAT, ge)

    entre_exp = 0
    worker_exp = 0
    MITupdate!(MIT_partial, AggrOutputE_new, AggrLabDemE_new, tax_rev_entre_new, g_tilde, g_tilde_e, SSall, A_t, Aentre_local_t, PartialPolW, PartialPolE, g, ge, TempW, TempE, P, G, Opt, entre_exp, worker_exp)

    MIT_partial.AggrOutputE_t .= AggrOutputE_new
    MIT_partial.AggrLabDemE_t .= AggrLabDemE_new
    MIT_partial.tax_rev_entre_t .= tax_rev_entre_new
    MIT_partial.C_t .= MIT_partial.ConsumpW_t .+ MIT_partial.ConsumpE_t
    MIT_partial.AggrLabDemC_t .= MIT_partial.AggrLabSup_t .- MIT_partial.AggrLabDemE_t
    MIT_partial.AggrIlliqSav_t[2:end] .= MIT_partial.AggrIlliqSavW_t[2:end] .+ MIT_partial.AggrIlliqSavE_t[2:end]
    MIT_partial.AggrLiqSav_t[2:end] .= MIT_partial.AggrLiqSavW_t[2:end] .+ MIT_partial.AggrLiqSavE_t[2:end]

    # build investment time series
    MIT_partial.AggrCapSup_t .= MIT_partial.AggrIlliqSav_t[1:(end-1)] .- MIT_partial.q_t
    MIT_partial.AggrCapDemC_t .= MIT_partial.AggrCapSup_t .- MIT_partial.AggrCapDemE_t
    MIT_partial.KW_t .= MIT_partial.AggrCapSup_t .* MIT_partial.AggrIlliqSavW_t[1:(end-1)] ./ MIT_partial.AggrIlliqSav_t[1:(end-1)]
    MIT_partial.KE_t .= MIT_partial.AggrCapSup_t .* MIT_partial.AggrIlliqSavE_t[1:(end-1)] ./ MIT_partial.AggrIlliqSav_t[1:(end-1)]
    # Keep the steady-state terminal boundary used by MIT_init.
    MIT_partial.InvC_M_t[end] = SSall.InvC_M
    MIT_partial.InvE_M_t[end] = SSall.InvE_M
    MIT_partial.Inv_M_t[end] = SSall.Inv_M
    MIT_partial.InvCW_t[end] = SSall.InvCW
    MIT_partial.InvCE_t[end] = SSall.InvCE
    for tt = (G.N_t-1):-1:1
        MIT_partial.InvC_M_t[tt] = (MIT_partial.AggrCapDemC_t[tt+1] - MIT_partial.AggrCapDemC_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT_partial.AggrCapDemC_t[tt]
        MIT_partial.InvE_M_t[tt] = (MIT_partial.AggrCapDemE_t[tt+1] - MIT_partial.AggrCapDemE_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT_partial.AggrCapDemE_t[tt]
        MIT_partial.Inv_M_t[tt] = MIT_partial.InvC_M_t[tt] .+ MIT_partial.InvE_M_t[tt]
        MIT_partial.InvCW_t[tt] = (MIT_partial.KW_t[tt+1] - MIT_partial.KW_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT_partial.KW_t[tt]
        MIT_partial.InvCE_t[tt] = (MIT_partial.KE_t[tt+1] - MIT_partial.KE_t[tt]) ./ G.dt[tt] .+ P.delta .* MIT_partial.KE_t[tt]
    end
    MIT_partial.InvE_t .= MIT_partial.InvE_M_t .+ MIT_partial.InvE_nonM_t
    MIT_partial.Inv_t .= MIT_partial.Inv_M_t .+ MIT_partial.InvE_nonM_t
    MIT_partial.OutputC_t .= MIT_partial.AggrCapDemC_t .^ P.alphaa .* MIT_partial.AggrLabDemC_t .^ (1.0 .- P.alphaa)
    MIT_partial.Output_total_t .= MIT_partial.OutputC_t .+ MIT_partial.AggrOutputE_t

    d_ela_MAT = reshape((PartialPolW[1].d .- PartialPolW[end].d), G.N_b, G.N_a, G.N_j, G.N_k) ./ G.a_grid_MAT_adj
    de_ela_MAT = reshape((PartialPolE[1].de .- PartialPolE[end].de), G.N_be, G.N_a, G.N_e, G.N_y) ./ G.a_grid_MAT_e_adj
    f_ela_MAT = reshape((PartialPolE[1].f .- PartialPolE[end].f), G.N_be, G.N_a, G.N_e, G.N_y) ./ G.e_grid_MAT_adj
    θ_ela_MAT = reshape((PartialPolE[1].θ .- PartialPolE[end].θ), G.N_be, G.N_a, G.N_e, G.N_y)

    # The following computes investment elasticities due to direct effect
    d_change_MAT = PartialPolW[1].d - PartialPolW[end].d
    de_change_MAT = PartialPolE[1].de - PartialPolE[end].de
    f_change_MAT = PartialPolE[1].f - PartialPolE[end].f
    θ_change_MAT = PartialPolE[1].θ - PartialPolE[end].θ

    println("Done")

    Partial = (MIT_partial=MIT_partial,
        PartialPolW=PartialPolW,
        PartialPolE=PartialPolE,
        f_ela_MAT=f_ela_MAT,
        d_ela_MAT=d_ela_MAT,
        de_ela_MAT=de_ela_MAT,
        θ_ela_MAT=θ_ela_MAT)

    if n_generators > 0
        entrepreneur_generators = [complete_entre_generator(Aentre_local_t[n], PartialPolE[n].θ, P, G, Opt) for n = 1:n_generators]
        generators = ( workers=A_t[1:n_generators], entrepreneurs=entrepreneur_generators,)
        return Partial, generators
    end

    return Partial

end

" `HHsolve_entre!(ve, PolE, Aentre, Bentre, RbE, rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants, adriftE, P, G, Opt, Delta_e, updatevalfunc, updatedecisions, TempE::EntreTemp)`

This file updates the value function using the HJB equation.
It performs one iteraton step using the implicit method and the upwind scheme.
Use `updatevalfunc` = false, to only update the transition matrix Aentre but not the value function.
Use `updatedecisions` = false, to only update the drifts in the transition matrix Aentre but not the decisions on d, f and c.
"
function HHsolve_entre!(ve, PolE, Aentre, Bentre, RbE, rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants, adriftE, P, G, Opt, Delta_e, updatevalfunc, updatedecisions, TempE::EntreTemp)

    ce_F = TempE.ce_F
    ce_B = TempE.ce_B
    θ_B = TempE.θ_B
    θ_F = TempE.θ_F
    θ_0 = TempE.θ_0
    θnew = TempE.θnew
    b = TempE.b
    ve_new = TempE.ve_new
    dBB = TempE.dBB
    dBF = TempE.dBF
    dFB = TempE.dFB
    dFF = TempE.dFF
    fBB = TempE.fBB
    fBF = TempE.fBF
    fFB = TempE.fFB
    fFF = TempE.fFF
    d_B = TempE.d_B
    d_F = TempE.d_F
    f_B = TempE.f_B
    f_F = TempE.f_F
    sf_B = TempE.sf_B
    sf_F = TempE.sf_F
    sce_B = TempE.sce_B
    sce_F = TempE.sce_F
    sde_B = TempE.sde_B
    sde_F = TempE.sde_F
    Ic_F = TempE.Ic_F
    Ic_B = TempE.Ic_B
    Ic_0 = TempE.Ic_0
    Id_F = TempE.Id_F
    Id_B = TempE.Id_B
    If_F = TempE.If_F
    If_B = TempE.If_B
    dB = TempE.dB
    dF = TempE.dF
    fB = TempE.fB
    fF = TempE.fF
    ve_b_F = TempE.ve_b_F
    ve_b_B = TempE.ve_b_B
    ve_a_F = TempE.ve_a_F
    ve_a_B = TempE.ve_a_B
    ve_e_F = TempE.ve_e_F
    ve_e_B = TempE.ve_e_B
    ve_e_F_raw = TempE.ve_e_F_raw
    ve_e_B_raw = TempE.ve_e_B_raw
    ve_ee = TempE.ve_ee
    edriftB = TempE.edriftB
    edriftF = TempE.edriftF
    edrift0 = TempE.edrift0
    a = (P.nu .* P.alphaaE) / (1 - P.nu * (1 - P.alphaaE))

    if updatedecisions

        @views @. ve_b_F[1:end-1, :, :, :] = max((ve[2:end, :, :, :] - ve[1:end-1, :, :, :]) / G.delta_b_MAT_e, Opt.v_b_min)
        @views @. ve_b_B[2:end, :, :, :] = ve_b_F[1:end-1, :, :, :]

        ve_b_B[1, :, :, :] .= max.(RbE[1] * G.b_grid_e[1] .+ TransE[1, :, :, :] .+ max.(P.intoliqshare .* (FirmIncomeConstants[1, :, :, :]
                                                           .+
                                                           (1 - P.tauE) .* RealGrossProfitsnoLevPerCap[1, :, :, :] .* G.e_grid_MAT[1, :, :, :] .* θ_0[1, :, :, :] .^ a
                                                           .-
                                                           (1.0 .- P.tauE) .* (θ_0[1, :, :, :] .- 1.0) .* G.e_grid_MAT[1, :, :, :] .* rk), 0.0), Opt.min_c) .^ (-P.sigma) # state constraint boundary condition
        ve_b_F[end, :, :, :] .= max.(RbE[end] * G.b_grid_e[end] .+ TransE[end, :, :, :] .+ max.(P.intoliqshare .* (FirmIncomeConstants[end, :, :, :] .+ (1 - P.tauE) .* RealGrossProfitsnoLevPerCap[end, :, :, :] .* G.e_grid_MAT[end, :, :, :] .* θ_0[end, :, :, :] .^ a .- (1.0 .- P.tauE) .* (θ_0[end, :, :, :] .- 1.0) .* G.e_grid_MAT[end, :, :, :] .* rk), 0.0), Opt.min_c) .^ (-P.sigma)   # no saving at upper end
        # Optimal consumption
        for xx = 1:G.Ne
            cF = ve_b_F[xx]^(-1.0 / P.sigma)
            cB = ve_b_B[xx]^(-1.0 / P.sigma)
            ce_F[xx] = max(cF, Opt.min_c)
            ce_B[xx] = max(cB, Opt.min_c)
        end

        @views @. ve_a_F[:, 1:end-1, :, :] = max((ve[:, 2:end, :, :] - ve[:, 1:end-1, :, :]) / G.delta_a_MAT_e, Opt.v_a_min)
        @views @. ve_a_B[:, 2:end, :, :] = ve_a_F[:, 1:end-1, :, :]

        @views @. ve_e_F_raw[:, :, 1:end-1, :] = (ve[:, :, 2:end, :] - ve[:, :, 1:end-1, :]) / G.delta_e_MAT
        @views @. ve_e_F[:, :, 1:end-1, :] = max(ve_e_F_raw[:, :, 1:end-1, :], Opt.v_e_min)
        @views @. ve_e_B_raw[:, :, 2:end, :] = ve_e_F_raw[:, :, 1:end-1, :]
        @views @. ve_e_B[:, :, 2:end, :] = max(ve_e_B_raw[:, :, 2:end, :], Opt.v_e_min)
        if Opt.leverage_mode == :risky
            # Store e*v_ee. Using the raw slopes is important: the positive
            # truncation used for investment choices would distort curvature.
            @views @. ve_ee[:, :, 2:end-1, :] =
                G.e_grid_MAT[:, :, 2:end-1, :] * 2.0 *
                (ve_e_F_raw[:, :, 2:end-1, :] - ve_e_B_raw[:, :, 2:end-1, :]) /
                (G.delta_e_MAT[:, :, 1:end-1, :] + G.delta_e_MAT[:, :, 2:end, :])

            # boundaries: simple one-sided approximations by copying neighbors
            ve_ee[:, :, 1, :] .= ve_ee[:, :, 2, :]
            ve_ee[:, :, G.N_e, :] .= ve_ee[:, :, G.N_e-1, :]

            damp = 0.9

            for yy = 1:G.N_y
                for aa = 1:G.N_a
                    for bb = 1:G.N_be
                        for ee = 2:G.N_e

                            if abs(lev_fun(θ_B[bb, aa, ee, yy], RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_B[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt)) < 1e-10
                                # do nothing, already solved

                                continue
                            elseif lev_fun(P.λ, RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_B[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt) >= 0.0
                                θnew[bb, aa, ee, yy] = P.λ

                            elseif lev_fun(1.0, RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_B[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt) <= 0.0
                                θnew[bb, aa, ee, yy] = 1.0

                            else
                                θnew[bb, aa, ee, yy] = bisection_intervalwidth(θ -> lev_fun(θ, RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_B[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt), 1.0, P.λ)
                            end

                            θ_B[bb, aa, ee, yy] = exp(log(θnew[bb, aa, ee, yy]) + damp * (log(θ_B[bb, aa, ee, yy]) - log(θnew[bb, aa, ee, yy]))) # damping in logs to ensure positivity
                        end
                    end
                end
            end
            for yy = 1:G.N_y
                for aa = 1:G.N_a
                    for bb = 1:G.N_be
                        for ee = 2:G.N_e
                            if abs(lev_fun(θ_F[bb, aa, ee, yy], RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_F[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt)) < 1e-10
                                # do nothing, already solved

                                continue
                            elseif lev_fun(P.λ, RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_F[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt) >= 0.0
                                θnew[bb, aa, ee, yy] = P.λ
                            elseif lev_fun(1.0, RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_F[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt) <= 0.0
                                θnew[bb, aa, ee, yy] = 1.0

                            else
                                θnew[bb, aa, ee, yy] = bisection_intervalwidth(θ -> lev_fun(θ, RealGrossProfitsnoLevPerCap[bb, aa, ee, yy], a, rk, P, ve_ee[bb, aa, ee, yy], ve_b_F[bb, aa, ee, yy], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt), 1.0, P.λ)
                            end

                            θ_F[bb, aa, ee, yy] = exp(log(θnew[bb, aa, ee, yy]) + damp * (log(θ_F[bb, aa, ee, yy]) - log(θnew[bb, aa, ee, yy]))) # damping in logs to ensure positivity
                        end
                    end
                end
            end
        elseif Opt.leverage_mode == :riskfree
            for xx in eachindex(θ_0)
                profit = RealGrossProfitsnoLevPerCap[xx]
                if profit > 0.0 && rk > 0.0
                    θstar = clamp((a * profit / rk)^(1.0 / (1.0 - a)),1.0,P.λ)
                elseif rk <= 0.0
                    θstar = P.λ
                else
                    θstar = 1.0
                end
                θ_B[xx] = θstar
                θ_F[xx] = θstar
                θ_0[xx] = θstar
            end
        end
        edriftB .= FirmIncomeConstants .+ (1 - P.tauE) .* RealGrossProfitsnoLevPerCap .* G.e_grid_MAT .* θ_B .^ a .- (1.0 .- P.tauE) .* (θ_B .- 1.0) .* G.e_grid_MAT .* rk
        edriftF .= FirmIncomeConstants .+ (1 - P.tauE) .* RealGrossProfitsnoLevPerCap .* G.e_grid_MAT .* θ_F .^ a .- (1.0 .- P.tauE) .* (θ_F .- 1.0) .* G.e_grid_MAT .* rk

        #Compute optimal deposits
        depodrift!(dBB, ve_a_B, ve_b_B, G.a_grid_MAT_e_adj, P.χ0, P.χ1, P.χ2)
        depodrift!(dFB, ve_a_B, ve_b_F, G.a_grid_MAT_e_adj, P.χ0, P.χ1, P.χ2)
        depodrift!(dBF, ve_a_F, ve_b_B, G.a_grid_MAT_e_adj, P.χ0, P.χ1, P.χ2)
        depodrift!(dFF, ve_a_F, ve_b_F, G.a_grid_MAT_e_adj, P.χ0, P.χ1, P.χ2)

        for xx = 1:G.Ne

            dBB[xx] = clamp(dBB[xx], -Opt.max_d, Opt.max_d)
            dBF[xx] = clamp(dBF[xx], -Opt.max_d, Opt.max_d)
            dFF[xx] = clamp(dFF[xx], -Opt.max_d, Opt.max_d)
            dFB[xx] = clamp(dFB[xx], -Opt.max_d, Opt.max_d)
        end
        #UPWIND SCHEME
        upwind!(d_B, dBF, dBB)
        #state constraints at amin and amax
        @views begin
            d_B[:, 1, :, :] .= (dBF[:, 1, :, :] .> 1e-12) .* dBF[:, 1, :, :] #make sure d>=0 at amin, don't use VaB(:,1,:)
            d_B[:, end, :, :] .= (dBB[:, end, :, :] .< -1e-12) .* dBB[:, end, :, :] #make sure d<=0 at amax, don't use VaF(:,J,:)
            d_B[1, 1, :, :] .= max.(d_B[1, 1, :, :], 0.0)
            d_B[end, :, :, :] .= max.(d_B[end, :, :, :], 0.0)
        end

        #Compute optimal firm investments
        depodrift!(fBB, ve_e_B, ve_b_B, G.e_grid_MAT_adj, 0.0, P.χe1, P.χe2)

        depodrift!(fFB, ve_e_B, ve_b_F, G.e_grid_MAT_adj, 0.0, P.χe1, P.χe2)

        depodrift!(fBF, ve_e_F, ve_b_B, G.e_grid_MAT_adj, 0.0, P.χe1, P.χe2)

        depodrift!(fFF, ve_e_F, ve_b_F, G.e_grid_MAT_adj, 0.0, P.χe1, P.χe2)

        for xx = 1:G.Ne

            fBB[xx] = clamp(fBB[xx], -Opt.max_d, Opt.max_d)
            fBF[xx] = clamp(fBF[xx], -Opt.max_d, Opt.max_d)
            fFF[xx] = clamp(fFF[xx], -Opt.max_d, Opt.max_d)
            fFB[xx] = clamp(fFB[xx], -Opt.max_d, Opt.max_d)
        end

        #UPWIND SCHEME
        upwind!(f_B, fBF, fBB)
        #state constraints at amin and amax
        @views begin
            f_B[:, :, 1, :] .= (fBF[:, :, 1, :] .> 1e-12) .* fBF[:, :, 1, :] #make sure d>=0 at amin, don't use VaB(:,1,:)
            f_B[:, :, end, :] .= (fBB[:, :, end, :] .< -1e-12) .* fBB[:, :, end, :] #make sure d<=0 at amax, don't use VaF(:,J,:)
            f_B[1, :, 1, :] .= max.(f_B[1, :, 1, :], 0.0)
            f_B[end, :, :, :] .= max.(f_B[end, :, :, :], 0.0)
        end
    end

    #split drift of b and upwind separately

    scebuild!(sce_B, RbE, TransE, ce_B, edriftB, FirmIncomeConstants, G, P, Opt)
    sdebuild!(sde_B, d_B, G.a_grid_MAT_e_adj, P.χ0, P.χ1, P.χ2)

    sdebuild!(sf_B, f_B, G.e_grid_MAT_adj, 0.0, P.χe1, P.χe2)

    if updatedecisions
        upwind!(d_F, dFF, dFB)

        @views begin
            d_F[end, :, :, :] .= max.(d_F[end, :, :, :], 0.0)
            #state constraints at amin and amax
            d_F[:, 1, :, :] .= (dFF[:, 1, :, :] .> 1e-12) .* dFF[:, 1, :, :] #make sure d>=0 at amin, don't use VaB(:,1,:)
            d_F[:, end, :, :] .= (dFB[:, end, :, :] .< -1e-12) .* dFB[:, end, :, :] #make sure d<=0 at amax, don't use VaF(:,J,:)
            f_F .= (fFF .> 0) .* fFF .+ (fFB .< 0.0) .* fFB
            f_F[end, :, :, :] .= max.(f_F[end, :, :, :], 0.0)
            #state constraints at amin and amax
            f_F[:, :, 1, :] .= (fFF[:, :, 1, :] .> 1e-12) .* fFF[:, :, 1, :] #make sure d>=0 at amin, don't use VaB(:,1,:)
            f_F[:, :, end, :] .= (fFB[:, :, end, :] .< -1e-12) .* fFB[:, :, end, :] #make sure d<=0 at amax, don't use VaF(:,J,:)
        end
    end
    #split drift of b and upwind separately

    sdebuild!(sde_F, d_F, G.a_grid_MAT_e_adj, P.χ0, P.χ1, P.χ2)

    scebuild!(sce_F, RbE, TransE, ce_F, edriftF, FirmIncomeConstants, G, P, Opt)
    sde_F[end, :, :, :] .= min.(sde_F[end, :, :, :], 0.0)

    sdebuild!(sf_F, f_F, G.e_grid_MAT_adj, 0.0, P.χe1, P.χe2)
    sf_F[end, :, :, :] .= min.(sf_F[end, :, :, :], 0)

    if updatedecisions
        driftgap2 = 1e-6
        Ic_B .= sce_B .< -driftgap2
        Ic_F .= (sce_F .> driftgap2) .& (.!Ic_B)
        Ic_0 .= (.!Ic_B) .& (.!Ic_F)

        Id_F .= sde_F .> driftgap2
        Id_B .= (sde_B .< -driftgap2) .& (.!Id_F)
        Id_B[1, :, :, :] .= 0
        Id_F[end, :, :, :] .= 0
        Id_B[end, :, :, :] .= 1 #don't use VbF at bmax so as not to pick up articial state constraint

        dB .= Id_B .* dBB .+ Id_F .* dFB
        dF .= Id_B .* dBF .+ Id_F .* dFF
        dB[:, 1, :, :] .= max.(dB[:, 1, :, :], 0.0)
        dF[end, :, :, :] .= max.(dF[end, :, :, :], 0.0)
        dB[end, :, :, :] .= max.(dB[end, :, :, :], 0.0)

        If_F .= sf_F .> driftgap2
        If_B .= (sf_B .< -driftgap2) .& (.!If_F)
        If_B[1, :, :, :] .= 0
        If_F[end, :, :, :] .= 0
        If_B[end, :, :, :] .= 1

        fB .= If_B .* fBB .+ If_F .* fFB
        fF .= If_B .* fBF .+ If_F .* fFF
        fB[:, :, 1, :] .= max.(fB[:, :, 1, :], 0.0)

        fF[end, :, :, :] .= max.(fF[end, :, :, :], 0.0)
        fB[end, :, :, :] .= max.(fB[end, :, :, :], 0.0)

        if Opt.leverage_mode == :risky
            ind = CartesianIndices(Ic_0)
            for xx = eachindex(Ic_0)
                bb = ind[xx][1]
                aa = ind[xx][2]
                ee = ind[xx][3]
                yy = ind[xx][4]
                if ee > 1 && Ic_0[xx]
                    if abs(lev_fun0(θ_0[xx], FirmIncomeConstants[xx], RealGrossProfitsnoLevPerCap[xx], a, rk, RbE[bb], G.b_grid_e[bb], G.e_grid[ee], TransE[xx], P, ve_ee[xx], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt)) < 1e-10
                        #do nothing, already solved

                        continue
                    elseif lev_fun0(P.λ, FirmIncomeConstants[xx], RealGrossProfitsnoLevPerCap[xx], a, rk, RbE[bb], G.b_grid_e[bb], G.e_grid[ee], TransE[xx], P, ve_ee[xx], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt) >= 0.0
                        θnew[bb, aa, ee, yy] = P.λ
                    elseif lev_fun0(1.0, FirmIncomeConstants[xx], RealGrossProfitsnoLevPerCap[xx], a, rk, RbE[bb], G.b_grid_e[bb], G.e_grid[ee], TransE[xx], P, ve_ee[xx], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt) <= 0.0
                        θnew[bb, aa, ee, yy] = 1.0
                    else
                        θnew[bb, aa, ee, yy] = bisection_intervalwidth(θ -> lev_fun0(θ, FirmIncomeConstants[xx], RealGrossProfitsnoLevPerCap[xx], a, rk, RbE[bb], G.b_grid_e[bb], G.e_grid[ee], TransE[xx], P, ve_ee[xx], Opt.leverage_mode, bb, aa, ee, yy, ve_e_F_raw, G, Opt), 1.0, P.λ)
                    end

                    θ_0[xx] = exp(log(θnew[bb, aa, ee, yy]) + damp * (log(θ_0[xx]) - log(θnew[bb, aa, ee, yy]))) # damping in logs to ensure positivity
                end

            end
        end
    end

    edrift0 .= FirmIncomeConstants .+ (1 - P.tauE) .* RealGrossProfitsnoLevPerCap .* G.e_grid_MAT .* θ_0 .^ a .- (1.0 .- P.tauE) .* (θ_0 .- 1.0) .* G.e_grid_MAT .* rk
    cebuild!(PolE.ce, RbE, Ic_0, G.b_grid_MAT_e, edrift0,
        FirmIncomeConstants, TransE, G, P, Opt)

    transit_mat_entre!(Aentre, Ic_B, Ic_F, Id_B, Id_F, If_B, If_F,
        sce_B, sce_F, sde_B, sde_F, sf_B, sf_F, dB, dF, fB, fF,
        adriftE, edrift0, FirmIncomeConstants, P, G, Opt, TempE)
    # For :none and :riskfree, diffusion and jump risk are already in G.Aexo.
    # For :risky, scale diffusion here and handle jumps matrix-free below.
    if Opt.leverage_mode == :risky

        for ii = eachindex(G.iadiff)
            i = G.iadiff[ii]
            pos = G.diff_nzidx[ii]
            Aentre.nzval[pos] += (θ_F[i] .^ 2.0 * Ic_F[i] + θ_B[i] .^ 2.0 * Ic_B[i] + θ_0[i] .^ 2.0 * Ic_0[i] - 1.0) .* G.vadiff[ii]
        end
    end
    bbuild!(b, G.Ne, PolE.ce, ce_B, ce_F, Ic_F, Ic_B, Ic_0, G.entre_lab_vec, ve, Delta_e, P) # this also sets up ce!
    jump_diagonal = 0.0

    if risky_jump_active(Opt, P, G)
        jump_diagonal = P.lambda_e
        for yy = 1:G.N_y
            for ee = 1:G.N_e
                e = G.e_grid[ee]
                for aa = 1:G.N_a
                    for bb = 1:G.N_be
                        θ = θ_F[bb, aa, ee, yy] * Ic_F[bb, aa, ee, yy] + θ_B[bb, aa, ee, yy] * Ic_B[bb, aa, ee, yy] + θ_0[bb, aa, ee, yy] * Ic_0[bb, aa, ee, yy]
                        expected_v = 0.0
                        for nn = eachindex(G.jump_u)
                            u = G.jump_u[nn]
                            enew = max(1.0 - θ * (1.0 - u), 0.0) * e
                            if enew <= G.e_grid[1]
                                vjump = ve[bb, aa, 1, yy]
                            elseif enew >= G.e_grid[end]
                                vjump = ve[bb, aa, end, yy]
                            else
                                pos = searchsortedfirst(G.e_grid, enew)
                                weightright = (enew - G.e_grid[pos-1]) / (G.e_grid[pos] - G.e_grid[pos-1])
                                weightleft = 1.0 - weightright
                                vjump = weightleft * ve[bb, aa, pos-1, yy] + weightright * ve[bb, aa, pos, yy]
                            end
                            expected_v += G.jump_weights[nn] * vjump
                        end
                        xx = G.linindices[bb, aa, ee, yy]
                        b[xx] += P.lambda_e * expected_v
                    end
                end
            end
        end
    end
    if updatevalfunc
        x = (1.0 / Delta_e + P.ρE + jump_diagonal)
        Bassemble!(Bentre, x, Aentre, G.iaexo, G.jaexo)
        linear_solve!(ve_new, Bentre, b, TempE,
            Opt.entrepreneur_HJB_solver, Opt.bicgstab_tol;
            label="entrepreneur HJB")
        crit = compcrit(ve_new, ve)
        ve .= reshape(ve_new, G.N_be, G.N_a, G.N_e, G.N_y)
    else
        crit = 0.0
    end

    PolE.θ .= vec(θ_F) .* vec(Ic_F) .+ vec(θ_B) .* vec(Ic_B) .+ vec(θ_0) .* vec(Ic_0)

    return crit

end

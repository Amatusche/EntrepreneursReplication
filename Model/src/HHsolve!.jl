function HHsolve!(v, Pol, A, B, Rb, w, Trans, lab0, adrift, P, G, Opt, Delta, Temp)

    ### This file updates the value function using the HJB equation.
    ### It performs one iteraton step using the implicit method and the upwind scheme

    c_0 = Temp.c_0
    c_F = Temp.c_F
    c_B = Temp.c_B
    lab_F = Temp.lab_F
    lab_B = Temp.lab_B
    lab = Temp.lab
    b = Temp.b
    v_new = Temp.v_new
    dBB = Temp.dBB
    dBF = Temp.dBF
    dFB = Temp.dFB
    dFF = Temp.dFF
    d_B = Temp.d_B
    d_F = Temp.d_F
    sc_B = Temp.sc_B
    sc_F = Temp.sc_F
    sd_B = Temp.sd_B
    sd_F = Temp.sd_F
    Ic_F = Temp.Ic_F
    Ic_B = Temp.Ic_B
    Ic_0 = Temp.Ic_0
    Id_F = Temp.Id_F
    Id_B = Temp.Id_B
    dB = Temp.dB
    dF = Temp.dF
    v_prime_b_F = Temp.v_prime_b_F
    v_prime_b_B = Temp.v_prime_b_B
    v_prime_a_F = Temp.v_prime_a_F
    v_prime_a_B = Temp.v_prime_a_B

    AA = 0.0
    BB = 0.0
    CC = 0.0

    @views @. v_prime_b_F[1:end-1, :, :, :] = max((v[2:end, :, :, :] - v[1:end-1, :, :, :]) / G.delta_b_MAT, Opt.v_b_min)
    @views @. v_prime_b_B[2:end, :, :, :] = v_prime_b_F[1:end-1, :, :, :]

    v_prime_b_B[1, :, :, :] .= max.((1.0 .- P.tau) .* w .* lab0[1, :, :, :] .* G.z_grid_MAT[1, :, :, :] .+ Rb[1] .* G.b_grid[1] .+ Trans[1, :, :, :], Opt.min_c) .^ (-P.sigma) # state constraint boundary condition
    v_prime_b_F[end, :, :, :] .= max.((1.0 .- P.tau) .* w .* lab0[end, :, :, :] .* G.z_grid_MAT[end, :, :, :] .+ Rb[end] .* G.b_grid[end] .+ Trans[end, :, :, :], Opt.min_c) .^ (-P.sigma)   # no saving at upper end

    v_prime_a_F .= 0.0
    v_prime_a_B .= 0.0

    @views @. v_prime_a_F[:, 1:end-1, :, :] = max((v[:, 2:end, :, :] - v[:, 1:end-1, :, :]) / G.delta_a_MAT, Opt.v_a_min)
    @views @. v_prime_a_B[:, 2:end, :, :] = v_prime_a_F[:, 1:end-1, :, :]

    for xx = 1:G.N
        if v_prime_b_F[xx] .> Opt.min_c
            c_F[xx] = v_prime_b_F[xx] .^ (-1.0 ./ P.sigma)
        else
            c_F[xx] = Opt.min_c .^ (-1.0 ./ P.sigma)
        end
        if v_prime_b_B[xx] .> Opt.min_c
            c_B[xx] = v_prime_b_B[xx] .^ (-1.0 ./ P.sigma)
        else
            c_B[xx] = Opt.min_c .^ (-1.0 ./ P.sigma)
        end
    end

    for xx = 1:G.N
        lab_F[xx] = ((1.0 .- P.tau) .* w .* G.z_grid_MAT[xx] ./ P.varphi .* v_prime_b_F[xx]) .^ P.gamma
        if lab_F[xx] > P.max_l
            lab_F[xx] = P.max_l
        end
        lab_B[xx] = ((1.0 .- P.tau) .* w .* G.z_grid_MAT[xx] ./ P.varphi .* v_prime_b_B[xx]) .^ P.gamma
        if lab_B[xx] > P.max_l
            lab_B[xx] = P.max_l
        end
    end

    #Compute optimal deposits
    depodrift!(dBB, v_prime_a_B, v_prime_b_B, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)
    depodrift!(dFB, v_prime_a_B, v_prime_b_F, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)
    depodrift!(dBF, v_prime_a_F, v_prime_b_B, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)
    depodrift!(dFF, v_prime_a_F, v_prime_b_F, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)

    for xx = 1:G.N
        if dBB[xx] > Opt.max_d
            dBB[xx] = Opt.max_d
            if Opt.display == 2
                println("             High optimal deposits")
            end
        end
        if dBF[xx] > Opt.max_d
            dBF[xx] = Opt.max_d
        end
        if dFF[xx] > Opt.max_d
            dFF[xx] = Opt.max_d
        end
        if dFB[xx] > Opt.max_d
            dFB[xx] = Opt.max_d
        end
    end
    #UPWIND SCHEME
    upwind!(d_B, dBF, dBB)
    #state constraints at amin and amax
    @views begin
        d_B[:, 1, :, :] .= (dBF[:, 1, :, :] .> 10.0^(-12.0)) .* dBF[:, 1, :, :] #make sure d>=0 at amin, don't use VaB(:,1,:)
        d_B[:, end, :, :] .= (dBB[:, end, :, :] .< -10.0^(-12.0)) .* dBB[:, end, :, :] #make sure d<=0 at amax, don't use VaF(:,J,:)
        d_B[1, 1, :, :] .= max.(d_B[1, 1, :, :], 0.0)
        d_B[end, :, :, :] .= max.(d_B[end, :, :, :], 0.0)
    end

    #split drift of b and upwind separately
    sc_B .= (1.0 .- P.tau) .* w .* G.z_grid_MAT .* lab_B .+ Rb .* G.b_grid_MAT .+ Trans .- c_B
    sdebuild!(sd_B, d_B, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)

    upwind!(d_F, dFF, dFB)
    @views begin
        d_F[end, :, :, :] .= max.(d_F[end, :, :, :], 0)
        #state constraints at amin and amax
        d_F[:, 1, :, :] .= (dFF[:, 1, :, :] .> 10^(-12)) .* dFF[:, 1, :, :] #make sure d>=0 at amin, don't use VaB(:,1,:)
        d_F[:, end, :, :] .= (dFB[:, end, :, :] .< -10^(-12)) .* dFB[:, end, :, :] #make sure d<=0 at amax, don't use VaF(:,J,:)
    end
    #split drift of b and upwind separately
    sc_F .= (1.0 .- P.tau) .* w .* G.z_grid_MAT .* lab_F .+ Rb .* G.b_grid_MAT .+ Trans .- c_F
    sdebuild!(sd_F, d_F, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)

    sd_F[end, :, :, :] .= min.(sd_F[end, :, :, :], 0)
    Ic_B .= sc_B .< (-10.0 .^ -12.0)
    Ic_F .= (sc_F .> (10.0 .^ -12.0)) .& (.!Ic_B)
    Ic_0 .= (.!Ic_B) .& (.!Ic_F)

    Id_F .= sd_F .> (10.0 .^ -12.0)
    Id_B .= (sd_B .< (-10.0 .^ -12.0)) .& (.!Id_F)
    Id_B[1, :, :, :] .= 0
    Id_F[end, :, :, :] .= 0
    Id_B[end, :, :, :] .= 1 #don't use VbF at bmax so as not to pick up articial state constraint

    dB .= Id_B .* dBB .+ Id_F .* dFB
    dF .= Id_B .* dBF .+ Id_F .* dFF
    dB[:, 1, :, :] .= max.(dB[:, 1, :, :], 0.0)
    dF[end, :, :, :] .= max.(dF[end, :, :, :], 0.0)
    dB[end, :, :, :] .= max.(dB[end, :, :, :], 0.0)

    x0 = 1.2
    ind = CartesianIndices(Ic_0)
    for xx = findall(vec(Ic_0))
        bb = ind[xx][1]
        if P.sigma .* P.gamma == 1
            BB = Rb[bb] .* G.b_grid[bb] .+ Trans[xx]
            AA = (1.0 .- P.tau) .* w .* G.z_grid_MAT[xx]
            CC = (AA ./ P.varphi) .^ (-P.sigma)
            lab0[xx] = -BB ./ (2.0 .* AA) .+ sqrt((BB ./ (2.0 .* AA)) .^ 2.0 .+ 1.0 ./ (CC .* AA))
            if AA .* lab0[xx] .+ BB < Opt.min_c
                lab0[xx] = 1.0 ./ (Opt.min_c .* CC)
            end
        else
            lab_eq(x) = lab_fun(x, P.tau, w, G.z_grid_MAT[xx], Rb[bb], G.b_grid[bb], Trans[xx], P.varphi, P.sigma, P.gamma, Opt.min_c)

            lab0[xx] = min.(fzero(lab_eq, x0), P.max_l)
        end
        c_0[xx] = max.((1.0 .- P.tau) .* w .* G.z_grid_MAT[xx] .* lab0[xx] .+ Rb[bb] .* G.b_grid_MAT[xx] .+ Trans[xx], Opt.min_c)
    end

    Ic_0 = (.!Ic_F) .& (.!Ic_B)
    for xx = 1:G.N
        if Ic_F[xx]
            Pol.c[xx] = c_F[xx]
        elseif Ic_B[xx]
            Pol.c[xx] = c_B[xx]
        else
            Pol.c[xx] = c_0[xx]
        end
    end
    for xx = 1:G.N
        lab[xx] = min.(((1.0 .- P.tau) .* w .* G.z_grid_MAT[xx] ./ P.varphi .* Pol.c[xx] .^ (-P.sigma)) .^ (P.gamma), P.max_l)
    end

    bbuild!(b, G.N, Pol.c, c_B, c_F, Ic_F, Ic_B, Ic_0, lab, v, Delta, P)

    transit_mat!(A, Ic_B, Ic_F, Id_B, Id_F, sd_B, sd_F, sc_B, sc_F, dB, dF, adrift, G, Temp)
    copyto!(v_new, v)
    x = (1.0 ./ Delta .+ P.ρ)
    Bassemble!(B, x, A, G.iaz, G.jaz)
    linear_solve!(v_new, B, b, Temp, Opt.worker_HJB_solver,
        Opt.bicgstab_tol; label="worker HJB")

    crit = compcrit(v_new, v)
    v .= reshape(v_new, G.N_b, G.N_a, G.N_j, G.N_k)

    Pol.lab .= vec(lab)

    return crit
end

function transit_mat!(A, Ic_B, Ic_F, Id_B, Id_F, sd_B, sd_F, sc_B, sc_F,
    dB, dF, adrift, G, Temp)

    Xb = Temp.trans_Xb
    Y = Temp.trans_Y
    Zb = Temp.trans_Zb
    Xa = Temp.trans_Xa
    Za = Temp.trans_Za
    backward_drift = Temp.trans_B
    forward_drift = Temp.trans_F

    # Liquid-asset transitions
    @views begin
        Xb .= -1.0 .* Ic_B[2:G.N] .* sc_B[2:G.N]
        Xb .-= Id_B[2:G.N] .* sd_B[2:G.N]
        Xb ./= G.rep_delta_b_B[2:G.N]

        Y .= (Ic_B[:] .* sc_B[:] .+ Id_B[:] .* sd_B[:]) ./ G.rep_delta_b_B .-
             (Ic_F[:] .* sc_F[:] .+ Id_F[:] .* sd_F[:]) ./ G.rep_delta_b_F

        Zb .= (Ic_F[1:G.N-1] .* sc_F[1:G.N-1] .+
               Id_F[1:G.N-1] .* sd_F[1:G.N-1]) ./ G.rep_delta_b_F[1:G.N-1]
    end

    # Passive illiquid-asset transitions
    backward_drift .= min.(dB, 0.0)
    forward_drift .= max.(dF, 0.0) .+ adrift
    @views forward_drift[:, end, :, :] .= 0.0

    @views begin
        Xa .= -1.0 .* backward_drift[G.N_b+1:G.N] ./ G.rep_delta_a_B[G.N_b+1:G.N]
        Y .+= backward_drift[:] ./ G.rep_delta_a_B .- forward_drift[:] ./ G.rep_delta_a_F
        Za .= forward_drift[1:G.N-G.N_b] ./ G.rep_delta_a_F[1:G.N-G.N_b]
    end

    Aassemble!(A, Y, Xa, Xb, Za, Zb, G)

    return nothing
end

"""
    transit_mat_fixed_worker!(A, Pol, Rb, w, Trans, adrift, P, G, Temp)

Build the worker transition matrix when consumption, labor, and deposit policies
are held fixed, but liquid-asset accumulation is evaluated at the realized
return `Rb`. This is used by `worker_exp == 3`: workers choose their policies
using the steady-state liquid return, while their balance sheets evolve at the
actual return.

Unlike `transit_mat!`, this routine obtains transition directions from the
realized drifts rather than reusing the upwind indicators from the perceived-
return problem.
"""
function transit_mat_fixed_worker!(A, Pol, Rb, w, Trans, adrift, P, G, Temp)
    Xb = Temp.trans_Xb
    Y = Temp.trans_Y
    Zb = Temp.trans_Zb
    Xa = Temp.trans_Xa
    Za = Temp.trans_Za

    consumption_drift = Temp.sc_B
    deposit_liquid_drift = Temp.sd_B
    backward_drift = Temp.trans_B
    forward_drift = Temp.trans_F

    c = reshape(Pol.c, G.N_b, G.N_a, G.N_j, G.N_k)
    lab = reshape(Pol.lab, G.N_b, G.N_a, G.N_j, G.N_k)
    d = reshape(Pol.d, G.N_b, G.N_a, G.N_j, G.N_k)

    # Realized liquid-account drifts. Policies remain those chosen under the
    # perceived return; only the budget flows use the actual return.
    @. consumption_drift = (1.0 - P.tau) * w * G.z_grid_MAT * lab +
        Rb * G.b_grid_MAT + Trans - c
    sdebuild!(deposit_liquid_drift, d, G.a_grid_MAT_adj, P.χ0, P.χ1, P.χ2)

    # Preserve the model's separate consumption/saving and deposit-flow
    # components, but determine both directions from their realized signs.
    @. backward_drift = min(consumption_drift, 0.0) +
        min(deposit_liquid_drift, 0.0)
    @. forward_drift = max(consumption_drift, 0.0) +
        max(deposit_liquid_drift, 0.0)
    @views begin
        backward_drift[1, :, :, :] .= 0.0
        forward_drift[end, :, :, :] .= 0.0

        Xb .= -backward_drift[2:G.N] ./ G.rep_delta_b_B[2:G.N]
        Y .= backward_drift[:] ./ G.rep_delta_b_B .-
             forward_drift[:] ./ G.rep_delta_b_F
        Zb .= forward_drift[1:G.N-1] ./ G.rep_delta_b_F[1:G.N-1]
    end

    # The deposit policy is fixed as well. Its realized movement in the
    # illiquid account is combined with the passive illiquid return.
    @. backward_drift = min(d, 0.0)
    @. forward_drift = max(d, 0.0) + adrift
    @views begin
        backward_drift[:, 1, :, :] .= 0.0
        forward_drift[:, end, :, :] .= 0.0

        Xa .= -backward_drift[G.N_b+1:G.N] ./ G.rep_delta_a_B[G.N_b+1:G.N]
        Y .+= backward_drift[:] ./ G.rep_delta_a_B .-
              forward_drift[:] ./ G.rep_delta_a_F
        Za .= forward_drift[1:G.N-G.N_b] ./ G.rep_delta_a_F[1:G.N-G.N_b]
    end

    Aassemble!(A, Y, Xa, Xb, Za, Zb, G)
    return nothing
end

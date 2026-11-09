"""Return whether policy-dependent entrepreneurial jump risk is active."""
function risky_jump_active(Opt, P, G)
    return Opt.leverage_mode == :risky &&
           Opt.entreriskjump &&
           P.lambda_e > 0.0 &&
           P.eta > 0.0 &&
           !isempty(G.jump_u)
end

"""
Derivative of the risky-jump HJB term with respect to leverage, divided by
current business equity. The derivative is taken through the same piecewise
linear interpolation used for the jump value in `HHsolve_entre!`.
"""
function risky_jump_foc(θ, bb, aa, ee, yy, ve_e_F_raw, P, G, Opt)
    risky_jump_active(Opt, P, G) || return 0.0

    e = G.e_grid[ee]
    e <= G.e_grid[1] && return 0.0

    derivative = 0.0
    for nn = eachindex(G.jump_u)
        destruction = 1.0 - G.jump_u[nn]
        remaining_share = max(1.0 - θ * destruction, 0.0)
        enew = remaining_share * e

        # The implemented jump map is constant once it reaches either grid
        # boundary, so its derivative is zero there.
        if enew <= G.e_grid[1] || enew >= G.e_grid[end]
            continue
        end

        pos = searchsortedfirst(G.e_grid, enew)
        value_slope = ve_e_F_raw[bb, aa, pos-1, yy]
        derivative -= G.jump_weights[nn] * destruction * value_slope
    end

    return P.lambda_e * derivative
end

"""
Construct the risky-leverage jump generator for a given leverage policy.

Rows are origin states. A jump with realization `u` maps equity to
`max(1 - theta * (1-u), 0) * e`. Post-jump mass is distributed linearly
between the two surrounding equity-grid points.
"""
function risky_jump_generator(theta, P, G)
    n_jump = length(G.jump_u)
    n_entries = G.Ne * (1 + 2 * n_jump)

    q_i = Int64[]
    q_j = Int64[]
    q_v = Float64[]
    sizehint!(q_i, n_entries)
    sizehint!(q_j, n_entries)
    sizehint!(q_v, n_entries)

    for yy = 1:G.N_y
        for ee = 1:G.N_e
            e = G.e_grid[ee]
            for aa = 1:G.N_a
                for bb = 1:G.N_be
                    i = G.linindices[bb, aa, ee, yy]

                    push!(q_i, i)
                    push!(q_j, i)
                    push!(q_v, -P.lambda_e)

                    for nn = eachindex(G.jump_u)
                        enew = max(1.0 - theta[i] * (1.0 - G.jump_u[nn]), 0.0) * e
                        quadrature_weight = G.jump_weights[nn]
                        if enew <= G.e_grid[1]
                            j = G.linindices[bb, aa, 1, yy]
                            push!(q_i, i)
                            push!(q_j, j)
                            push!(q_v, P.lambda_e * quadrature_weight)
                        elseif enew >= G.e_grid[end]
                            j = G.linindices[bb, aa, G.N_e, yy]
                            push!(q_i, i)
                            push!(q_j, j)
                            push!(q_v, P.lambda_e * quadrature_weight)
                        else
                            pos = searchsortedfirst(G.e_grid, enew)
                            weightright = (enew - G.e_grid[pos-1]) / (G.e_grid[pos] - G.e_grid[pos-1])
                            weightleft = 1.0 - weightright
                            jleft = G.linindices[bb, aa, pos-1, yy]
                            jright = G.linindices[bb, aa, pos, yy]
                            push!(q_i, i)
                            push!(q_j, jleft)
                            push!(q_v, P.lambda_e * quadrature_weight * weightleft)
                            push!(q_i, i)
                            push!(q_j, jright)
                            push!(q_v, P.lambda_e * quadrature_weight * weightright)
                        end
                    end
                end
            end
        end
    end

    return sparse(q_i, q_j, q_v, G.Ne, G.Ne)
end

"""Return a complete entrepreneur generator, including risky jumps if active."""
function complete_entre_generator(Aentre_local, theta, P, G, Opt)
    if risky_jump_active(Opt, P, G)
        return Aentre_local + risky_jump_generator(theta, P, G)
    end
    return copy(Aentre_local)
end

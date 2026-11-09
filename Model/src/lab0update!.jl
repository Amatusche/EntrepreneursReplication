function lab0update!(lab0, Trans, w, Rb, G, P, Opt)
    AA = 0.0
    BB = 0.0
    CC = 0.0
    x0 = 0.0
    for bb = [1 G.N_b]
        for jj = 1:G.N_j
            for kk = 1:G.N_k
                if P.sigma .* P.gamma == 1
                    BB = Rb[bb] .* G.b_grid[bb] .+ Trans[bb, 1, jj, kk]
                    AA = (1.0 .- P.tau) .* w .* G.z_grid[jj, kk]
                    CC = (AA ./ P.varphi) .^ (-P.sigma)
                    lab0[bb, :, jj, kk] .= -BB ./ (2.0 .* AA) .+ sqrt.((BB ./ (2.0 .* AA)) .^ 2.0 .+ 1.0 ./ (CC .* AA))
                    if (AA .* lab0[bb, 1, jj, kk] .+ BB) < Opt.min_c
                        lab0[bb, :, jj, kk] .= 1.0 ./ (Opt.min_c .* CC)
                    end
                else
                    x0 = lab0[bb, 1, jj, kk]
                    lab_eq(x) = lab_fun(x, P.tau, w, G.z_grid[jj, kk], Rb[bb], G.b_grid[bb], Trans[bb, 1, jj, kk], P.varphi, P.sigma, P.gamma, Opt.min_c)
                    lab0[bb, :, jj, kk] .= min.(fzero(lab_eq, x0), P.max_l)
                end
            end
        end
    end
end

"""
    fractional_grid_index(grid, value)

Return the fractional one-based index that linearly interpolates `value`
between two adjacent points of a strictly increasing `grid`.
"""
function fractional_grid_index(grid::AbstractVector{<:Real}, value::Real)
    isempty(grid) && throw(ArgumentError("grid must not be empty"))
    all(diff(grid) .> 0.0) ||
        throw(ArgumentError("grid must be strictly increasing"))

    if value < first(grid) || value > last(grid)
        throw(ArgumentError(
            "requested value $value is outside grid range " *
            "[$(first(grid)), $(last(grid))]"))
    end

    right_index = searchsortedfirst(grid, value)
    if right_index == 1 || grid[right_index] == value
        return Float64(right_index)
    end

    left_index = right_index - 1
    right_weight = (value - grid[left_index]) /
                   (grid[right_index] - grid[left_index])
    return left_index + right_weight
end

function grids(P, Opt; N_e=45)
    ## Import and set gridsizes
    N_j = 3 # number of productivity component z1 states
    N_k = 11 # number of productivity component z2 states
    N_y = 2
    N_a = 40
    e_curve = 0.25
    a_curve = 0.15
    N_t = 200 # number of transition periods
    @set! P.inheritance_index_e = max(min(P.inheritance_index_e, N_e), 1)
    a_max = 2000.0
    e_max = 2000.0
    N_b_pos = 40
    N_b_neg = 10
    N_b = N_b_pos + N_b_neg
    b_max = 1000.0
    N_be_pos = 40
    N_be_neg = 10
    N_be = N_be_pos + N_be_neg
    be_max = 1000.0
    ejumpgridpoints = 10 # Number of points used to discretize the continuous capital-destruction shock
    N_ab = N_a * N_b # number of asset states for workers
    N = N_ab * N_j * N_k # total number of states for workers
    N_abe = N_a * N_be # number of asset states for entrepreneurs
    N_bae = N_a * N_be * N_e # number of asset states for entrepreneurs
    Ne = N_bae * N_y # total number of states for entrepreneurs
    linindices = LinearIndices((N_be, N_a, N_e, N_y))

    ## I.a Liquid asset grid (b) for workers
    # Use same as in HANK
    b_grid = Vector{Float64}(undef, N_b)

    b_grid[(end-N_b_pos+1):end] = get_dt(N_b_pos, 0.28, 0, b_max)
    b_grid[1:(div(N_b_neg, 2)+1)] = get_dt(div(N_b_neg, 2) + 1, 0.4, P.underline_b, (P.underline_b + b_grid[N_b_neg+1]) / 2)
    for bb = (div(N_b_neg, 2)+2):N_b_neg
        b_grid[bb] = b_grid[N_b_neg+1] - (b_grid[N_b_neg+2-bb] - b_grid[1])
    end

    b_grid_MAT = repeat(b_grid, 1, N_a, N_j, N_k)
    delta_b_F = b_grid[2:end] - b_grid[1:(end-1)]
    delta_b_F = [delta_b_F; delta_b_F[end]]
    delta_b_B = b_grid[2:end] - b_grid[1:(end-1)]
    delta_b_B = [delta_b_B[1]; delta_b_B]
    delta_b_MAT = repeat(delta_b_F[1:(N_b-1)], 1, N_a, N_j, N_k)
    db_tilde = 0.5 * (delta_b_B + delta_b_F)
    db_tilde[1] = 0.5 * delta_b_F[1]
    db_tilde[end] = 0.5 * delta_b_B[end]
    b_low_pos = N_b_neg + 1   # find the index of the 0 on the b grid
    b_low_pos_e = N_be_neg + 1   # find the index of the 0 on the b grid for entrepreneurs

    ## I.b Liquid asset grid (b) for entrepreneurs
    if N_be == 1
        b_grid_e = [0.0]
        b_grid_MAT_e = repeat(b_grid_e, 1, N_a, N_e, N_y)
        delta_be_F = [1.0]
        delta_be_B = [1.0]
        delta_b_MAT_e = repeat(delta_be_F[1:(N_be-1)], 1, N_a, N_e, N_y)
        db_tilde_e = [1.0]
        be_low_pos = 1
    else
        b_grid_e = Vector{Float64}(undef, N_be)

        b_grid_e[(end-N_be_pos+1):end] = get_dt(N_be_pos, 0.28, 0, be_max)
        b_grid_e[1:(div(N_be_neg, 2)+1)] = get_dt(div(N_be_neg, 2) + 1, 0.4, P.underline_be, (P.underline_be + b_grid_e[N_be_neg+1]) / 2)
        for bb = (div(N_be_neg, 2)+2):N_be_neg
            b_grid_e[bb] = b_grid_e[N_be_neg+1] - (b_grid_e[N_be_neg+2-bb] - b_grid_e[1])
        end
        b_grid_MAT_e = repeat(b_grid_e, 1, N_a, N_e, N_y)
        delta_be_F = b_grid_e[2:end] - b_grid_e[1:(end-1)]
        delta_be_F = [delta_be_F; delta_be_F[end]]
        delta_be_B = b_grid_e[2:end] - b_grid_e[1:(end-1)]
        delta_be_B = [delta_be_B[1]; delta_be_B]
        delta_b_MAT_e = repeat(delta_be_F[1:(N_be-1)], 1, N_a, N_e, N_y)
        db_tilde_e = 0.5 * (delta_be_B + delta_be_F)
        db_tilde_e[1] = 0.5 * delta_be_F[1]
        db_tilde_e[end] = 0.5 * delta_be_B[end]
        be_low_pos = N_be_neg + 1   # find the index of the 0 on the b grid
    end

    ## II. Labor productivity grid (z): Poisson Jump-Drift Process
    # Use same as in HANK
    lambda_1 = 7.987252454850691E-002
    lambda_2 = 6.560658609206301E-003
    var_z1 = 1.73522438989472^2 # variance of new z1 draws
    var_z2 = 1.52585125038250^2 # variance of new z2 draws
    beta_1 = 0.761606393760105
    beta_2 = 9.391785115626633E-003

    # Import HANK discretization (z grids, phi and asset grids)
    HANK_ljump1 = readdlm("../in/ljump1.txt")
    HANK_ljump2 = readdlm("../in/ljump2.txt")
    phi_1 = HANK_ljump1 / lambda_1 + I(N_j)
    phi_1 = phi_1[1, :]
    phi_2 = HANK_ljump2 / lambda_2 + I(N_k)
    phi_2 = phi_2[1, :]
    z_1 = readdlm("../in/y1grid.txt")'
    z_2 = readdlm("../in/y2grid.txt")'
    adjfrischfrac = 0.85
    z_1_adj = z_1 / (1 + adjfrischfrac * P.gamma)
    z_2_adj = z_2 / (1 + adjfrischfrac * P.gamma)

    if N_j == 1
        delta_z1_F = [1.0]
        delta_z1_B = [1.0]
        z_1_stat = [1.0]'
    else
        delta_z1_F = z_1_adj[2:end] .- z_1_adj[1:(end-1)]
        delta_z1_F = [delta_z1_F; delta_z1_F[end]]
        delta_z1_B = z_1_adj[2:end] .- z_1_adj[1:(end-1)]
        delta_z1_B = [delta_z1_B[1]; delta_z1_B]
        z_1_stat = find_z_stat(N_j, beta_1, z_1_adj, phi_1, lambda_1, delta_z1_B, delta_z1_F)
    end
    if N_k == 1
        delta_z2_F = [1.0]
        delta_z2_B = [1.0]
        z_2_stat = [1.0]'
    else
        delta_z2_F = z_2_adj[2:end] .- z_2_adj[1:(end-1)]
        delta_z2_F = [delta_z2_F; delta_z2_F[end]]
        delta_z2_B = z_2_adj[2:end] .- z_2_adj[1:(end-1)]
        delta_z2_B = [delta_z2_B[1]; delta_z2_B]
        z_2_stat = find_z_stat(N_k, beta_2, z_2_adj, phi_2, lambda_2, delta_z2_B, delta_z2_F)
    end

    z_stat = z_1_stat .* z_2_stat'

    AvgLabProduct = 0.997047141890806
    exp_z_1 = exp.(z_1_adj) ./ (exp.(z_1_adj) * z_1_stat)
    exp_z_2 = exp.(z_2_adj) ./ (exp.(z_2_adj) * z_2_stat)
    z_grid = AvgLabProduct * exp_z_1' * exp_z_2

    z_grid_MAT = reshape(z_grid, 1, 1, N_j, N_k) .* ones(N_b, N_a, 1, 1)
    Az = exostate_transit_mat(N_ab, N_j, N_k, beta_1, beta_2, lambda_1, lambda_2, phi_1, phi_2, z_1_adj, z_2_adj, delta_z1_F, delta_z1_B, delta_z2_F, delta_z2_B)

    Apattern = spdiagm(-N_b => zeros((N_b+1):N), -1 => zeros(2:N),
        0 => zeros(N), 1 => zeros(1:(N-1)), N_b => zeros(1:(N-N_b)))
    i1, j1, vaz1 = findnz(Apattern)
    i2, j2, vaz2 = findnz(Az)
    iaz = [i1; i2]
    jaz = [j1; j2]
    vaz = [vaz1; vaz2]
    Az = sparse(iaz, jaz, vaz)
    iaz, jaz, vaz = findnz(Az)

    # Construct Lambda matrix capturing movement along z_1 and z_2
    # dimensions. This is used for "explicit part" of separated procedure.
    if Opt.sep_inc_states == 1
        Lambda = Matrix(undef, N_j, N_k)
        for jj = 1:N_j
            for kk = 1:N_k
                jk_part = (N_ab*N_j*(kk-1)+N_ab*(jj-1)+1):(N_ab*N_j*(kk-1)+N_ab*(jj-1)+N_ab)
                Lambda[jj, kk] = Az[jk_part, :]
            end
        end
    end

    b_vec = zeros(N)
    for jj = 1:N_j
        for kk = 1:N_k
            b_vec[(jj .- 1) .* N_b .* N_a .+ (kk .- 1) .* N_b .* N_a .* N_j .+ b_low_pos] = z_stat[jj, kk]
        end
    end
    b_vec .= b_vec ./ sum(b_vec)

    ## III. Illiquid asset grid (a)
    if N_a > 1
        a_grid = get_dt(N_a, a_curve, 0.0, a_max)
        for aa = 1:9
            a_grid[aa] = (aa - 1) * a_grid[10] / (10 - 1)  # procedure as for a_grid in HANK
        end
    else
        a_grid = [0.0]
    end

    if N_a == 1
        delta_a_F = [1.0]
        delta_a_B = [1.0]
    else
        delta_a_F = a_grid[2:end] - a_grid[1:(end-1)]
        delta_a_F = [delta_a_F; delta_a_F[end]]
        delta_a_B = a_grid[2:end] - a_grid[1:(end-1)]
        delta_a_B = [delta_a_B[1]; delta_a_B]
        delta_a_F = max.(delta_a_F, 10e-5)
        delta_a_B = max.(delta_a_B, 10e-5)
    end
    da_tilde = 0.5 * (delta_a_B + delta_a_F)
    da_tilde[1] = 0.5 * delta_a_F[1]
    da_tilde[end] = 0.5 * delta_a_B[end]
    a_grid_MAT = permutedims(repeat(a_grid, 1, N_b, N_j, N_k), (2, 1, 3, 4))
    delta_a_MAT = permutedims(repeat(delta_a_F[1:(N_a-1)], 1, N_b, N_j, N_k), [2, 1, 3, 4])
    a_grid_MAT_e = permutedims(repeat(a_grid, 1, N_be, N_e, N_y), [2, 1, 3, 4])
    delta_a_MAT_e = permutedims(repeat(delta_a_F[1:(N_a-1)], 1, N_be, N_e, N_y), [2, 1, 3, 4])
    a_grid_MAT_e_adj = max.(a_grid_MAT_e, Opt.a_min2)
    a_grid_MAT_adj = max.(a_grid_MAT, Opt.a_min2)

    ## IV. Business capital grid (e)
    jump_u = Float64[]
    jump_weights = Float64[]
    if N_e > 1
        e_grid = get_dt(N_e, e_curve, 0, e_max)
        nopoints = 0
        for ee = 1:nopoints
            e_grid[ee] = (ee - 1) * e_grid[nopoints+1] / nopoints  # procedure as for a_grid in HANK
        end

        e_grid_MAT = permutedims(repeat(e_grid, 1, N_be, N_a, N_y), [2, 3, 1, 4])
        delta_e_F = e_grid[2:end] - e_grid[1:(end-1)]
        delta_e_F = [delta_e_F; delta_e_F[end]]
        delta_e_B = e_grid[2:end] - e_grid[1:(end-1)]
        delta_e_B = [delta_e_B[1]; delta_e_B]
        delta_e_F = max.(delta_e_F, 10e-5)
        delta_e_B = max.(delta_e_B, 10e-5)
        delta_e_MAT = permutedims(repeat(delta_e_F[1:(N_e-1)], 1, N_be, N_a, N_y), [2, 3, 1, 4])
        de_tilde = 0.5 * (delta_e_B + delta_e_F)
        de_tilde[1] = 0.5 * delta_e_F[1]
        de_tilde[end] = 0.5 * delta_e_B[end]
        de_stacked = repeat(kron(de_tilde, ones(N_be * N_a, 1)), N_y, 1)[:, 1]
        grid_diag_e = spdiagm(0 => de_stacked)

        e_grid_MAT_adj = max.(e_grid_MAT, Opt.a_min2)

        ## V. Build transition matrix that captures entrepreneurial capital quality shock
        Adiff = spzeros(Ne, Ne)
        Ajump0 = spzeros(Ne, Ne)
        if Opt.entreriskjump && P.eta > 0 && P.lambda_e > 0
            # Capital-destruction jumps arrive as a Poisson process with rate lambda_e.
            # Conditional on a jump, e' = e*u, where u = exp(-epsilon) and
            # epsilon has an exponential distribution with rate eta (scale 1/eta).

            # E[u] = eta/(eta+1), so the mean fraction destroyed is 1/(eta+1).
            # Form equal-probability bins of epsilon = -log(u).
            pL = (0:(ejumpgridpoints-1)) ./ ejumpgridpoints
            a = quantile.(Exponential(1 / P.eta), pL)
            pU = (1:ejumpgridpoints) ./ ejumpgridpoints
            b = quantile.(Exponential(1 / P.eta), pU)
            # Represent each bin by E[u | a <= epsilon < b].
            # This conditional mean equals eta/(eta+1) times
            # (exp(-(eta+1)a) - exp(-(eta+1)b)) /
            # (exp(-eta*a) - exp(-eta*b)).
            jump_u = (P.eta / (P.eta + 1.0)) .* (exp.(-(P.eta + 1) .* a) .- exp.(-(P.eta + 1) .* b)) ./ (exp.(-P.eta .* a) .- exp.(-P.eta .* b))
            jump_weights = fill(1.0 / ejumpgridpoints, ejumpgridpoints)
            q_i = Int64[]
            q_j = Int64[]
            q_v = Float64[]
            for ee = 1:N_e
                push!(q_i, ee)
                push!(q_j, ee)
                push!(q_v, -P.lambda_e)

                for nn = 1:ejumpgridpoints
                    enew = jump_u[nn] * e_grid[ee] # new grid points after shock, for each original grid point
                    # find the position of enew on the grid
                    if enew <= e_grid[1]
                        push!(q_i, ee)
                        push!(q_j, 1)
                        push!(q_v, P.lambda_e * jump_weights[nn])
                    else
                        pos = searchsortedfirst(e_grid, enew)

                        weightright = (enew - e_grid[pos-1]) / (e_grid[pos] - e_grid[pos-1])
                        weightleft = 1 - weightright
                        push!(q_i, ee)
                        push!(q_j, pos)
                        push!(q_v, P.lambda_e * jump_weights[nn]* weightright)
                        push!(q_i, ee)
                        push!(q_j, pos-1)
                        push!(q_v, P.lambda_e * jump_weights[nn]* weightleft)
                    end
                end
            end
            Qe = sparse(q_i, q_j, q_v, N_e, N_e)
            # Lift to full entrepreneur state (ordering has e-blocks of size N_abe within each y block)
            Ajump0 = kron(sparse(I, N_y, N_y), kron(Qe, sparse(I, N_abe, N_abe)))
        end
        if Opt.entreriskdiffusion

            xi = e_grid .^ 2 * P.sigma_e^2 ./ (delta_e_F .* (delta_e_F + delta_e_B))
            nu = -e_grid .^ 2 * P.sigma_e^2 ./ (delta_e_F .* delta_e_B)
            zeta = e_grid .^ 2 * P.sigma_e^2 ./ (delta_e_B .* (delta_e_F + delta_e_B))
            updiag = zeta[1] * ones(N_a * N_be)
            for ee = 2:(N_e-1)
                updiag = [updiag; zeta[ee] * ones(N_a * N_be)]
            end
            #This will be the center diagonal of ExoTrans
            centdiag = (nu[1] + xi[1]) * ones(N_a * N_be)
            for ee = 2:(N_e-1)
                centdiag = [centdiag; nu[ee] * ones(N_a * N_be)]
            end
            centdiag = [centdiag; (nu[N_e] + zeta[N_e]) * ones(N_a * N_be)]
            #This will be the lower diagonal of ExoTrans
            lodiag = (xi[2]) * ones(N_a * N_be)
            for ee = 3:N_e
                lodiag = [lodiag; xi[ee] * ones(N_a * N_be)]
            end
            #Add up the upper, center, and lower diagonal into a sparse matrix
            for yy = 1:N_y

                Adiff[((yy-1)*N_e*N_a*N_be+1):(yy*N_e*N_a*N_be), ((yy-1)*N_e*N_a*N_be+1):(yy*N_e*N_a*N_be)] = spdiagm(0 => centdiag, -N_a * N_be => lodiag, N_a * N_be => updiag)
            end
        end
    else
        e_grid = [e_max]
        e_grid_MAT = permutedims(repeat(e_grid, 1, N_be, N_a, N_y), [2, 3, 1, 4])
        delta_e_F = [1.0]
        delta_e_F = [1.0]
        delta_e_B = [1.0]
        delta_e_B = [1.0]
        delta_e_F = [1.0]
        delta_e_B = [1.0]
        delta_e_MAT = ones(N_be, N_a, N_e, N_y)
        de_tilde = 0.5 * (delta_e_B + delta_e_F)
        de_tilde[1] = 0.5 * delta_e_F[1]
        de_tilde[end] = 0.5 * delta_e_B[end]
        de_stacked = repeat(kron(de_tilde, ones(N_be * N_a, 1)), N_y, 1)[:, 1]
        grid_diag_e = spdiagm(0 => de_stacked)
        e_grid_MAT_adj = max.(e_grid_MAT, Opt.a_min2)
        Adiff = spzeros(Ne, Ne)
        Ajump0 = spzeros(Ne, Ne)
    end

    ## VI. Entrepreneurial talent
    if N_y > 1
        y_stat = find_y_stat(P.lambda_y_lh, P.lambda_y_hl)
        f(x) = get_y_states(x, y_stat, P.log_y_spread)[1]

        y_low = Optim.optimize(f, -5.0, 5.0).minimizer
        y_high = get_y_states(y_low, y_stat, P.log_y_spread)[2]

        E_log_y = y_stat * [y_low y_high]'
        Var_log_y = y_stat * ([y_low y_high]' .- E_log_y) .^ 2
        y_grid = exp.([y_low y_high])

    else
        y_stat = [1]
        y_grid = [1]
    end

    Ay = exostate_transit_mat_entre(N_bae, N_y, P.lambda_y_lh, P.lambda_y_hl)

    Aepattern = spdiagm(-N_be * N_a => zeros((N_be*N_a+1):Ne), -N_be => zeros((N_be+1):Ne), -1 => zeros(2:Ne),
        0 => zeros(Ne), 1 => zeros(1:(Ne-1)), N_be => zeros(1:(Ne-N_be)), N_be * N_a => zeros(1:(Ne-N_be*N_a)))
    i1, j1, vaexo1 = findnz(Aepattern)
    Abase = Ay + Adiff
    if Opt.leverage_mode != :risky
        Abase += Ajump0
    end
    i2, j2, vaexo2 = findnz(Abase)
    iaexo = [i1; i2]
    jaexo = [j1; j2]
    vaexo = [vaexo1; vaexo2]
    Aexo = sparse(iaexo, jaexo, vaexo)
    iaexo, jaexo, vaexo = findnz(Aexo)
    iadiff, jadiff, vadiff = findnz(Adiff)
    diff_nzidx = sparse_nz_positions(Aexo, iadiff, jadiff)

    b_vecE = zeros(Ne)

    inheritance_index_be = N_be_neg + 1
    inheritance_index_a = 1
    inheritance_index_e_left = floor(Int, P.inheritance_index_e)
    inheritance_index_e_right = ceil(Int, P.inheritance_index_e)
    right_weight = P.inheritance_index_e - inheritance_index_e_left
    left_weight = 1.0 - right_weight
    inherited_business_equity =
        left_weight * e_grid[inheritance_index_e_left] +
        right_weight * e_grid[inheritance_index_e_right]
    inh_wealth = P.share_entre .* P.frac_inh_wealth * (b_grid_e[inheritance_index_be] + a_grid[inheritance_index_a] +
                                                       inherited_business_equity)
    for yy = 1:N_y
        b_vecE[(yy-1)*N_bae+0*(N_be*N_a)+b_low_pos_e] = y_stat[yy] * (1 - P.frac_inh_wealth)
        b_vecE[(yy-1)*N_bae+(inheritance_index_e_left-1)*(N_be*N_a)+(inheritance_index_a-1)*N_be+inheritance_index_be] =
            b_vecE[(yy-1)*N_bae+(inheritance_index_e_left-1)*(N_be*N_a)+(inheritance_index_a-1)*N_be+inheritance_index_be] +
            y_stat[yy] * P.frac_inh_wealth * left_weight
        b_vecE[(yy-1)*N_bae+(inheritance_index_e_right-1)*(N_be*N_a)+(inheritance_index_a-1)*N_be+inheritance_index_be] =
            b_vecE[(yy-1)*N_bae+(inheritance_index_e_right-1)*(N_be*N_a)+(inheritance_index_a-1)*N_be+inheritance_index_be] +
            y_stat[yy] * P.frac_inh_wealth * right_weight
    end

    rep_deltaE_b_B = repeat(delta_be_B, N_a * N_e * N_y)
    rep_deltaE_b_F = repeat(delta_be_F, N_a * N_e * N_y)
    rep_deltaE_a_B = repeat(kron(delta_a_B, ones(N_be)), N_e * N_y)
    rep_deltaE_a_F = repeat(kron(delta_a_F, ones(N_be)), N_e * N_y)
    rep_deltaE_e_B = repeat(kron(delta_e_B, ones(N_be * N_a)), N_y)
    rep_deltaE_e_F = repeat(kron(delta_e_F, ones(N_be * N_a)), N_y)

    rep_delta_b_B = repeat(delta_b_B, N_a * N_j * N_k)
    rep_delta_b_F = repeat(delta_b_F, N_a * N_j * N_k)
    rep_delta_a_B = repeat(kron(delta_a_B, ones(N_b)), N_j * N_k)
    rep_delta_a_F = repeat(kron(delta_a_F, ones(N_b)), N_j * N_k)

    ## MIT shock
    # Construct time grid here (in HANK they use N_t = 200 and P.deltatransmax = 40. N_t = 100 and P.deltatransmax = 80 achieves roughly the same time horizon (ca 2000
    # quarters))
    deltatransparam = 0.35 # all as in HANK, parameters for time grid
    deltatransmin = 1 / 3
    deltatransmax = 10
    dt = get_dt(N_t, deltatransparam, deltatransmin, deltatransmax)   # time per subperiod (fraction of a quarter we are looking at in each subperiod), power-spaced grid

    shock_ac = exp(-0.5)
    mp_shock_t = shock(-0.01, shock_ac, N_t, dt)

    entre_lab_vec = P.entre_lab * ones(Ne)

    G = (N_j=N_j,
        N_y=N_y,
        N_a=N_a,
        N_e=N_e,
        N_k=N_k,
        a_max=a_max,
        e_max=e_max,
        N_b=N_b,
        N_ab=N_ab,
        N_abe=N_abe,
        N_be=N_be,
        N=N,
        N_bae=N_bae,
        Ne=Ne,
        b_grid=b_grid,
        b_grid_MAT=b_grid_MAT,
        b_grid_e=b_grid_e,
        b_grid_MAT_e=b_grid_MAT_e,
        delta_b_F=delta_b_F,
        delta_b_B=delta_b_B,
        delta_b_MAT=delta_b_MAT,
        delta_b_MAT_e=delta_b_MAT_e,
        db_tilde=db_tilde,
        b_low_pos=b_low_pos, lambda_1=lambda_1,
        lambda_2=lambda_2,
        beta_1=beta_1,
        beta_2=beta_2,
        z_stat=z_stat, AvgLabProduct=AvgLabProduct,
        z_grid=z_grid, z_grid_MAT=z_grid_MAT,
        Az=Az,
        iaz=iaz,
        jaz=jaz,
        vaz=vaz, b_vec=b_vec,
        a_grid=a_grid,
        delta_a_F=delta_a_F,
        delta_a_B=delta_a_B,
        da_tilde=da_tilde,
        a_grid_MAT=a_grid_MAT,
        delta_a_MAT=delta_a_MAT,
        a_grid_MAT_e=a_grid_MAT_e,
        delta_a_MAT_e=delta_a_MAT_e,
        a_grid_MAT_e_adj=a_grid_MAT_e_adj,
        a_grid_MAT_adj=a_grid_MAT_adj, e_grid=e_grid,
        e_grid_MAT=e_grid_MAT,
        delta_e_F=delta_e_F,
        delta_e_B=delta_e_B,
        delta_e_MAT=delta_e_MAT,
        de_tilde=de_tilde, e_grid_MAT_adj=e_grid_MAT_adj, rep_deltaE_b_B=rep_deltaE_b_B,
        rep_deltaE_b_F=rep_deltaE_b_F,
        rep_deltaE_a_B=rep_deltaE_a_B,
        rep_deltaE_a_F=rep_deltaE_a_F,
        rep_deltaE_e_B=rep_deltaE_e_B,
        rep_deltaE_e_F=rep_deltaE_e_F, rep_delta_b_B=rep_delta_b_B,
        rep_delta_b_F=rep_delta_b_F,
        rep_delta_a_B=rep_delta_a_B,
        rep_delta_a_F=rep_delta_a_F, y_grid=y_grid,
        y_stat=y_stat,
        Ay=Ay,
        Adiff=Adiff,
        Aexo=Aexo,
        iaexo=iaexo,
        jaexo=jaexo,
        vaexo=vaexo,
        iadiff=iadiff,
        jadiff=jadiff,
        vadiff=vadiff,
        diff_nzidx=diff_nzidx,
        jump_weights=jump_weights,
        jump_u=jump_u, b_vecE=b_vecE, entre_lab_vec=entre_lab_vec,

        inheritance_index_be=inheritance_index_be,
        inheritance_index_a=inheritance_index_a,
        inh_wealth=inh_wealth,
        dt=dt,
        dt_cum=cumsum(dt),
        mp_shock_t=mp_shock_t,
        N_t=N_t,
        linindices=linindices)

end

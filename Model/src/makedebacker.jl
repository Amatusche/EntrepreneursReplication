"""
    debacker_mobility_stats(transition_matrix; min_decile_loss=3)

Compute mobility statistics from a square, row-stochastic transition matrix whose
rows and columns are ordered income groups. A large loss is a transition for which
`origin - destination >= min_decile_loss`. A large gain is defined symmetrically
as `destination - origin >= min_decile_loss`.

The statistics average over origin groups with equal weight, as appropriate for a
decile transition matrix. `rank_autocorrelation` is the correlation between origin
and destination group ranks under the joint distribution implied by the matrix.
`avg_deciles_lost` is conditional on a large loss.
"""
function debacker_mobility_stats(transition_matrix; min_decile_loss=3)
    n_origin, n_destination = size(transition_matrix)
    n_origin == n_destination || throw(ArgumentError("transition matrix must be square"))
    min_decile_loss >= 1 || throw(ArgumentError("min_decile_loss must be positive"))

    row_sums = sum(transition_matrix; dims=2)
    all(isapprox.(row_sums, 1.0; atol=1e-10)) ||
        throw(ArgumentError("transition-matrix rows must sum to one"))

    large_loss_probability = 0.0
    large_gain_probability = 0.0
    weighted_deciles_lost = 0.0
    mean_origin = (n_origin + 1) / 2
    mean_destination = sum(destination * transition_matrix[origin, destination] / n_origin
        for origin in axes(transition_matrix, 1), destination in axes(transition_matrix, 2))
    rank_covariance = 0.0
    origin_rank_variance = 0.0
    destination_rank_variance = 0.0
    for origin in axes(transition_matrix, 1)
        for destination in axes(transition_matrix, 2)
            probability = transition_matrix[origin, destination] / n_origin
            origin_deviation = origin - mean_origin
            destination_deviation = destination - mean_destination
            rank_covariance += probability * origin_deviation * destination_deviation
            origin_rank_variance += probability * origin_deviation^2
            destination_rank_variance += probability * destination_deviation^2

            deciles_lost = origin - destination
            if deciles_lost >= min_decile_loss
                large_loss_probability += probability
                weighted_deciles_lost += deciles_lost * probability
            end
            if -deciles_lost >= min_decile_loss
                large_gain_probability += probability
            end
        end
    end

    large_loss_probability > 0.0 ||
        throw(ArgumentError("transition matrix contains no losses of at least $min_decile_loss groups"))
    origin_rank_variance > 0.0 && destination_rank_variance > 0.0 ||
        throw(ArgumentError("transition matrix implies zero rank variance"))
    avg_deciles_lost = weighted_deciles_lost / large_loss_probability
    rank_autocorrelation = rank_covariance /
        sqrt(origin_rank_variance * destination_rank_variance)

    return (
        immobility_ratio=mean(diag(transition_matrix)),
        rank_autocorrelation=rank_autocorrelation,
        large_loss_probability=large_loss_probability,
        large_gain_probability=large_gain_probability,
        avg_deciles_lost=avg_deciles_lost,
    )
end

"""
    makedebacker(includeexits)

Set up the published DeBacker, Panousi, and Ramnath (2023) one-year business-
income transition matrix. If `includeexits` is true, business exits are assigned to
the first decile; otherwise exits are removed and each row is renormalized.
"""
function makedebacker(includeexits)

    DeBackerDataFull = [
        0.857 0.017 0.022 0.022 0.023 0.019 0.016 0.011 0.006 0.004 0.003;
        0.122 0.476 0.146 0.052 0.036 0.031 0.023 0.023 0.026 0.031 0.033;
        0.163 0.139 0.332 0.150 0.067 0.046 0.036 0.021 0.019 0.019 0.008;
        0.235 0.050 0.132 0.303 0.121 0.063 0.035 0.025 0.017 0.015 0.005;
        0.217 0.033 0.059 0.137 0.251 0.139 0.070 0.042 0.029 0.016 0.007;
        0.159 0.030 0.045 0.065 0.147 0.254 0.147 0.074 0.049 0.022 0.009;
        0.132 0.024 0.035 0.039 0.070 0.149 0.278 0.159 0.072 0.034 0.009;
        0.086 0.022 0.023 0.025 0.045 0.070 0.152 0.327 0.176 0.060 0.015;
        0.050 0.026 0.017 0.021 0.025 0.039 0.064 0.166 0.378 0.182 0.033;
        0.033 0.030 0.014 0.013 0.014 0.020 0.028 0.052 0.160 0.489 0.147;
        0.015 0.035 0.007 0.005 0.006 0.007 0.008 0.014 0.028 0.123 0.752]

    if includeexits
        # Only exclude inflows from zero state (no business income)
        DeBackerDataFull = DeBackerDataFull[2:end, :]
    else
        # Exclude zero state (no business income) alltogether
        DeBackerDataFull = DeBackerDataFull[2:end, :]
        DeBackerDataFull = DeBackerDataFull[:, 2:end]
    end

    for nn = eachindex(DeBackerDataFull[:,1])
        DeBackerDataFull[nn, :] = DeBackerDataFull[nn, :] / sum(DeBackerDataFull[nn, :])
    end

    if includeexits
        # Add exits and transitions into first decile
        DeBackerDataFull[:, 2] = DeBackerDataFull[:, 2] + DeBackerDataFull[:, 1]
        # Remove first column with exits
        DeBackerDataFull = DeBackerDataFull[:, 2:end]
    end

    # Transform to quintiles
    DeBackerDataQuint = Array{Float64}(undef, 5, 5)
    n = 1
    m = 1
    for nn = 1:5
        for mm = 1:5
            DeBackerDataQuint[nn, mm] = sum(DeBackerDataFull[n:n+1, m:m+1])
            m = m + 2
        end
        m = 1
        n = n + 2
    end
    # Normalize
    for nn = 1:size(DeBackerDataQuint, 2)
        DeBackerDataQuint[nn, :] = DeBackerDataQuint[nn, :] / sum(DeBackerDataQuint[nn, :])
    end

    return DeBackerDataQuint, DeBackerDataFull
end

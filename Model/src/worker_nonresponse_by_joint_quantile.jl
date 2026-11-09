"""
    worker_nonresponse_by_joint_quantile(g, ge, networth_w, networth_e,
                                         share_entre; n_groups=10)

Split workers into responding and nonresponding groups using quantiles of the
joint household net-worth distribution. Within each quantile, select one
nonresponding worker for every entrepreneur, so the aggregate mass and the
net-worth-bin distribution of nonresponding workers match those of
entrepreneurs. Probability mass at a discrete net-worth cutoff is split across
adjacent quantiles.
"""
function worker_nonresponse_by_joint_quantile(
    g, ge, networth_w, networth_e, share_entre; n_groups::Integer=10)

    n_groups > 0 || throw(ArgumentError("n_groups must be positive"))
    0.0 <= share_entre < 1.0 ||
        throw(ArgumentError("share_entre must lie in [0, 1)"))
    length(g) == length(networth_w) ||
        throw(DimensionMismatch("worker distribution and net worth differ in length"))
    length(ge) == length(networth_e) ||
        throw(DimensionMismatch("entrepreneur distribution and net worth differ in length"))

    minimum(g) >= -1e-12 ||
        throw(ArgumentError("worker distribution contains negative mass"))
    minimum(ge) >= -1e-12 ||
        throw(ArgumentError("entrepreneur distribution contains negative mass"))

    g_normalized = max.(g, 0.0)
    ge_normalized = max.(ge, 0.0)
    worker_total = sum(g_normalized)
    entrepreneur_total = sum(ge_normalized)
    worker_total > 0.0 || error("worker distribution has zero mass")
    entrepreneur_total > 0.0 || error("entrepreneur distribution has zero mass")
    g_normalized ./= worker_total
    ge_normalized ./= entrepreneur_total

    worker_mass_by_nw = Dict{Float64, Float64}()
    entrepreneur_mass_by_nw = Dict{Float64, Float64}()
    for ii in eachindex(g_normalized, networth_w)
        nw = Float64(networth_w[ii])
        mass = (1.0 - share_entre) * g_normalized[ii]
        if mass > 0.0
            worker_mass_by_nw[nw] = get(worker_mass_by_nw, nw, 0.0) + mass
        end
    end
    for ii in eachindex(ge_normalized, networth_e)
        nw = Float64(networth_e[ii])
        mass = share_entre * ge_normalized[ii]
        if mass > 0.0
            entrepreneur_mass_by_nw[nw] =
                get(entrepreneur_mass_by_nw, nw, 0.0) + mass
        end
    end

    nw_values = sort!(collect(union(keys(worker_mass_by_nw),
                                    keys(entrepreneur_mass_by_nw))))
    total_mass = sum(values(worker_mass_by_nw)) +
                 sum(values(entrepreneur_mass_by_nw))
    worker_mass_by_group = zeros(n_groups)
    entrepreneur_mass_by_group = zeros(n_groups)

    cumulative_mass = 0.0
    for nw in nw_values
        worker_mass = get(worker_mass_by_nw, nw, 0.0)
        entrepreneur_mass = get(entrepreneur_mass_by_nw, nw, 0.0)
        point_mass = worker_mass + entrepreneur_mass
        next_cumulative_mass = cumulative_mass + point_mass
        for group in 1:n_groups
            lower = total_mass * (group - 1) / n_groups
            upper = total_mass * group / n_groups
            overlap = max(0.0, min(next_cumulative_mass, upper) -
                               max(cumulative_mass, lower))
            if overlap > 0.0
                fraction = overlap / point_mass
                worker_mass_by_group[group] += fraction * worker_mass
                entrepreneur_mass_by_group[group] +=
                    fraction * entrepreneur_mass
            end
        end
        cumulative_mass = next_cumulative_mass
    end

    group_mass = worker_mass_by_group .+ entrepreneur_mass_by_group
    all(isapprox.(group_mass, total_mass / n_groups;
                  rtol=1e-11, atol=1e-12)) ||
        error("joint net-worth groups do not have equal mass")
    nonresponse_fraction_by_group = similar(group_mass)
    for group in eachindex(group_mass)
        worker_mass_by_group[group] > 0.0 ||
            error("joint net-worth group $group contains no workers")
        nonresponse_fraction_by_group[group] =
            entrepreneur_mass_by_group[group] / worker_mass_by_group[group]
        nonresponse_fraction_by_group[group] <= 1.0 + 1e-12 ||
            error("joint net-worth group $group contains fewer workers than " *
                  "entrepreneurs, so equal-mass matching is infeasible")
        nonresponse_fraction_by_group[group] =
            min(nonresponse_fraction_by_group[group], 1.0)
    end
    matched_worker_mass_by_group = worker_mass_by_group .*
                                   nonresponse_fraction_by_group
    all(isapprox.(matched_worker_mass_by_group, entrepreneur_mass_by_group;
                  rtol=1e-11, atol=1e-12)) ||
        error("worker and entrepreneur masses do not match within each group")

    nonresponse_fraction_by_nw = Dict{Float64, Float64}()
    cumulative_mass = 0.0
    for nw in nw_values
        point_mass = get(worker_mass_by_nw, nw, 0.0) +
                     get(entrepreneur_mass_by_nw, nw, 0.0)
        next_cumulative_mass = cumulative_mass + point_mass
        nonresponse_fraction = 0.0
        for group in 1:n_groups
            lower = total_mass * (group - 1) / n_groups
            upper = total_mass * group / n_groups
            overlap = max(0.0, min(next_cumulative_mass, upper) -
                               max(cumulative_mass, lower))
            nonresponse_fraction += overlap / point_mass *
                                    nonresponse_fraction_by_group[group]
        end
        nonresponse_fraction_by_nw[nw] =
            clamp(nonresponse_fraction, 0.0, 1.0)
        cumulative_mass = next_cumulative_mass
    end

    gstartNA_raw = similar(g_normalized)
    for ii in eachindex(g_normalized, networth_w)
        gstartNA_raw[ii] = g_normalized[ii] *
                           get(nonresponse_fraction_by_nw,
                               Float64(networth_w[ii]), 0.0)
    end
    shareNA = sum(gstartNA_raw)
    0.0 < shareNA < 1.0 ||
        error("nonresponding worker share must lie strictly between zero and one")

    gstartNA = gstartNA_raw ./ shareNA
    gstart = (g_normalized .- gstartNA_raw) ./ (1.0 - shareNA)
    minimum(gstart) >= -1e-12 ||
        error("responding worker distribution contains negative mass")
    gstart .= max.(gstart, 0.0)
    gstart ./= sum(gstart)

    reconstructed_g = (1.0 - shareNA) .* gstart .+ shareNA .* gstartNA
    reconstruction_error = maximum(abs.(reconstructed_g .- g_normalized))
    reconstruction_error <= 1e-12 ||
        error("worker split does not reconstruct the steady-state distribution")

    nonresponding_worker_mass = (1.0 - share_entre) * shareNA
    entrepreneur_mass = share_entre
    isapprox(nonresponding_worker_mass, entrepreneur_mass;
             rtol=1e-11, atol=1e-12) ||
        error("nonresponding worker and entrepreneur masses do not match")

    diagnostics = (; worker_mass_by_group, entrepreneur_mass_by_group,
                   matched_worker_mass_by_group,
                   nonresponse_fraction_by_group, group_mass,
                   nonresponding_worker_mass, entrepreneur_mass,
                   reconstruction_error)
    return gstart, gstartNA, shareNA, diagnostics
end

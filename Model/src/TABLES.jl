"Create the model tables used by the current paper draft."

using Printf, JLD2

include("HANKModel.jl")
using .HANKModel

const CURRENT_OPT = options()
const MODEL_PATH = CURRENT_OPT.outpath
const TABLE_PATH = MODEL_PATH * "Tables/"
mkpath(TABLE_PATH)

SSall, SSStats, P, G = load(joinpath(MODEL_PATH, "SS.jld2"),
    "SSall", "SSStats", "P", "G")
Targets = HANKModel.calibration_targets()
# Read model moments from this steady state, including its simulated mobility
# matrix, so tables do not depend on files left over from an optimizer run.
mobility = debacker_mobility_stats(SSStats.DeBackerModelFull)
CalibStats = merge((; (name => getfield(SSStats, name) for name in fieldnames(typeof(SSStats)))...),
    (EInc_immobility_ratio=mobility.immobility_ratio,
     EInc_large_loss_probability=mobility.large_loss_probability,
     EInc_large_gain_probability=mobility.large_gain_probability,
     EInc_avg_deciles_lost=mobility.avg_deciles_lost))
MIT = load(joinpath(MODEL_PATH, "MIT_Base.jld2"), "MIT")
MIT.converged[1] || error("MIT_Base.jld2 has not converged")
length(MIT.Output_total_t) == G.N_t ||
    error("MIT_Base.jld2 uses a different time grid than SS.jld2")

write_table(filename, contents) = write(TABLE_PATH * filename, contents)
fmt(x, digits=2) = @sprintf("%.*f", digits, x)

# Internally calibrated parameters, ordered and grouped as in the paper.
# Report transition intensities and pure preference discount rates;
# stored discount rates include mortality, and calibration uses their ratio.
io = IOBuffer()
println(io, "\\begin{tabular}{lcll}")
println(io, "\\toprule")
println(io, "Parameter & Value & Description & Main identifying target \\\\")
println(io, "\\midrule")
println(io, "\\multicolumn{4}{l}{Entrepreneurial income risk} \\\\")
println(io, "\$\\sigma_k\$ & ", fmt(P.sigma_e, 3), " & capital depreciation risk & income immobility ratio \\\\")
println(io, "\$y_h/y_l\$ & ", fmt(exp(P.log_y_spread), 3), " & productivity ratio & avg. deciles lost given large loss \\\\")
println(io, "\$\\lambda_{lh}\$ & ", fmt(P.lambda_y_lh, 3), " & low-to-high transition rate & probability of large income gain \\\\")
println(io, "\$\\lambda_{hl}\$ & ", fmt(P.lambda_y_hl, 3), " & high-to-low transition rate & probability of large income loss \\\\")
println(io, "\\midrule")
println(io, "\\multicolumn{4}{l}{Aggregate targets} \\\\")
println(io, "\$\\rho^W\$ & ", fmt(P.ρ - P.zeta, 4), " & worker discount rate & liquid assets/output \\\\")
println(io, "\$\\kappa\$ & ", fmt(P.kappa, 4), " & borrowing wedge & illiquid assets/output \\\\")
println(io, "\$\\chi^a_1\$ & ", fmt(P.χ1, 3), " & portfolio adjustment cost & wealthy hand-to-mouth share \\\\")
println(io, "\$\\chi^e_1\$ & ", fmt(P.χe1, 3), " & business-capital adjustment cost & entrepreneur htm share \\\\")
println(io, "\$\\bar y\$ & ", fmt(P.entre_prod, 3), " & average business productivity & private-business employment \\\\")
println(io, "\$\\rho\$ & ", fmt(P.ρE - P.zeta, 4), " & entrepreneur discount rate & relative average net worth \\\\")
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
if CURRENT_OPT.leverage_mode == :riskfree
    write_table("Table_C.1_CalibrationInternal.tex", String(take!(io)))
else
    write_table("Table_3_CalibrationInternal.tex", String(take!(io)))
end

# Externally calibrated parameters.
io = IOBuffer()
println(io, "\\begin{tabular}{lccll}")
println(io, "\\toprule")
println(io, "\\multicolumn{2}{c}{Parameter} & Value & Description & Source or target \\\\")
println(io, "\\midrule")
println(io, "Demographics & \$s_e\$ & ", fmt(100P.share_entre, 1), "\\% & entrepreneur share & SCF \\\\")
println(io, " & \$\\zeta\$ & ", fmt(P.zeta, 4), " & death rate & average lifetime 45 years \\\\")
println(io, "\\midrule")
println(io, "Preferences & \$\\varphi\$ & ", fmt(P.varphi, 1), " & labor disutility & average working time 8h/day \\\\")
println(io, " & \$\\gamma\$ & ", fmt(P.gamma, 1), " & inverse Frisch elasticity & KMV \\\\")
println(io, "\\midrule")
println(io, "Technology & \$\\underline b\$ & ", fmt(abs(P.underline_b), 1), " & borrowing limit & average quarterly labor income \\\\")
println(io, " & \$\\alpha\$ & ", fmt(P.alphaa, 2), " & capital share & KMV \\\\")
println(io, " & \$\\omega\$ & ", fmt(P.omega_liq, 2), " & dividend ratio & equal to capital share \\\\")
println(io, " & \$\\epsilon\$ & ", fmt(P.epsilon, 0), " & demand elasticity & KMV \\\\")
println(io, " & \$\\theta\$ & ", fmt(P.theta, 0), " & price adjustment cost & KMV \\\\")
println(io, " & \$\\delta\$ & ", fmt(P.delta, 3), " & depreciation rate & 7\\% p.a. \\\\")
println(io, " & \$\\nu\$ & ", fmt(P.nu, 2), " & returns to scale & \\textcite{Tan2021} \\\\")
println(io, " & \$\\chi^a_2\$ & ", fmt(P.χ2, 3), " & portfolio adjustment curvature & \\textcite{Alves2020} \\\\")
println(io, "\\midrule")
println(io, "Policy & \$\\phi\$ & ", fmt(P.phi, 2), " & inflation response & KMV \\\\")
if P.tauE == P.tau
    println(io, " & \$\\tau \$ & ", fmt(P.tau, 1), " & labor- and business-income tax & KMV \\\\")
else
    println(io, " & \$\\tau_l\$ & ", fmt(P.tau, 1), " & labor-income tax rate & KMV \\\\")
    println(io, " & \$\\tau_e\$ & ", fmt(P.tauE, 1), " & business-income tax rate & same as labor income \\\\")
end
println(io, " & \$T_t\$ & \$0.06Y_t\$ & lump-sum transfers & KMV \\\\")
println(io, " & \$\\bar r\$ & ", fmt(SSall.rb, 4), " & steady-state liquid rate & 2\\% p.a. \\\\")
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
write_table("Table_2_CalibrationExternal.tex", String(take!(io)))

# Current calibration targets and model fit, grouped as in the paper.
# Each row specifies its display multiplier and decimal precision explicitly.
target_groups = [
    ("Entrepreneurial income risk", [
        ("Income immobility ratio (\\%)", :EInc_immobility_ratio, 100.0, 1),
        ("Probability of large income loss (\\%)", :EInc_large_loss_probability, 100.0, 1),
        ("Probability of large income gain (\\%)", :EInc_large_gain_probability, 100.0, 1),
        ("Average deciles lost given a large loss", :EInc_avg_deciles_lost, 1.0, 2),
    ]),
    ("Aggregate targets", [
        ("Liquid assets / annual GDP (\$B/Y\$)", :AggrLiqSavAnnGDP, 1.0, 2),
        ("Illiquid assets / annual GDP (\$A/Y\$)", :AggrIlliqSavAnnGDP, 1.0, 2),
        ("Wealthy hand-to-mouth households (\\%)", :wHtmShare, 100.0, 1),
        ("Hand-to-mouth entrepreneurs (\\%)", :HtmShareE, 100.0, 1),
        ("Private-business employment share (\\%)", :LaborShareE, 100.0, 1),
        ("Average net worth: entrepreneurs / workers", :Avg_NW_ratio, 1.0, 2),
    ]),
]
io = IOBuffer()
println(io, "\\begin{tabular}{lrr}")
println(io, "\\toprule")
println(io, "Moment & Target & Model \\\\")
for (group, rows) in target_groups
    println(io, "\\midrule")
    println(io, "\\multicolumn{3}{l}{", group, "} \\\\")
    for (label, field, multiplier, digits) in rows
        println(io, label, " & ", fmt(multiplier * getfield(Targets, field), digits),
            " & ", fmt(multiplier * getfield(CalibStats, field), digits), " \\\\")
    end
end
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
if CURRENT_OPT.leverage_mode == :riskfree
    write_table("Table_C.2_TargetsVsModel.tex", String(take!(io)))
else
    write_table("Table_4_TargetsVsModel.tex", String(take!(io)))
end

# Untargeted net-worth moments.
io = IOBuffer()
println(io, "\\begin{tabular}{ccccc}")
println(io, "\\toprule")
println(io, " & Bottom 50\\% & Top 20\\% & Top 10\\% & Top 1\\% \\\\")
println(io, "\\midrule")
println(io, "Model & ", fmt(100(1 - SSStats.NWShares[1, 2]), 1), " & ",
    fmt(100SSStats.NWShares[2, 2], 1), " & ", fmt(100SSStats.NWShares[3, 2], 1),
    " & ", fmt(100SSStats.NWShares[6, 2], 1), " \\\\")
println(io, "SCF 2019 & 0.1 & 87.4 & 76.5 & 37.2 \\\\")
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
write_table("Table_6_WealthIneqVsData.tex", String(take!(io)))

entrepreneur_rows = [1, 2, 3, 4, 5, 6, 9]
io = IOBuffer()
println(io, "\\begin{tabular}{cccccccc}")
println(io, "\\toprule")
println(io, " & All & 1.+2. Q. & 3. Q. & 4. Q. & 5. Q. & Top 10\\% & Top 1\\% \\\\")
println(io, "\\midrule")
println(io, "Model & ", join((fmt(100SSStats.ShareEByNWGroup[i, 3], 1) for i in entrepreneur_rows), " & "), " \\\\")
println(io, "SCF 2019 & 7.1 & 1.5 & 4.2 & 8.9 & 19.4 & 25.6 & 48.6 \\\\")
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
write_table("Table_5_EntreShareVsData.tex", String(take!(io)))

# Decile transition matrices used to assess entrepreneurial income risk.
# Data exclude the zero-income state and are row-normalized, as in calibration.
# Differences are model minus data in percentage points, with a symmetric
# red/blue shading scale. Only the already-required xcolor package is needed.
function debacker_decile_table(matrix; differences=false)
    size(matrix) == (10, 10) || error("Expected a 10-by-10 decile matrix")
    all(isfinite, matrix) || error("Non-finite DeBacker transition entries")
    max_difference = differences ? maximum(abs, 100 .* matrix) : 0.0
    io = IOBuffer()
    println(io, "% Rows: initial income decile; columns: income decile one year later.")
    println(io, differences ? "% Model minus data, in percentage points." :
        "% Empirical transition probabilities, in percent; zero-income state excluded.")
    println(io, "\\begingroup")
    println(io, "\\setlength{\\tabcolsep}{3pt}")
    if differences
        println(io, "\\setlength{\\fboxsep}{1.5pt}")
        println(io, "\\definecolor{debackerpositive}{RGB}{178,24,43}")
        println(io, "\\definecolor{debackernegative}{RGB}{33,102,172}")
        println(io, "% Red: overprediction; blue: underprediction. Maximum absolute difference: ",
            fmt(max_difference, 3), " percentage points.")
    end
    println(io, "\\begin{tabular}{r*{10}{r}}")
    println(io, "\\toprule")
    println(io, "From / To & ", join(1:10, " & "), " \\\\")
    println(io, "\\midrule")
    for row in 1:10
        cells = String[]
        for col in 1:10
            value = 100 * matrix[row, col]
            if differences
                # Display rounded zeros neutrally, avoiding misleading +/-0.0.
                value = round(value; digits=1)
                label = value == 0 ? "0.0" : @sprintf("%+.1f", value)
                intensity = max_difference == 0 ? 0 :
                    round(Int, 65 * min(abs(value) / max_difference, 1.0))
                color = value > 0 ? "debackerpositive" : "debackernegative"
                fill = value == 0 ? "white" : color * "!" * string(intensity) * "!white"
                push!(cells, "\\colorbox{" * fill * "}{\\makebox[2.4em][r]{\\strut\$" * label * "\$}}")
            else
                push!(cells, fmt(value, 1))
            end
        end
        println(io, row, " & ", join(cells, " & "), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    println(io, "\\endgroup")
    return String(take!(io))
end

_, DeBackerDataFull = makedebacker(false)
DeBackerModelFull = SSStats.DeBackerModelFull
size(DeBackerModelFull) == size(DeBackerDataFull) == (10, 10) ||
    error("Model and empirical DeBacker matrices must both contain ten deciles")
for (name, matrix) in (("Data", DeBackerDataFull), ("Model", DeBackerModelFull))
    all(x -> x >= 0, matrix) || error("$name transition probabilities must be nonnegative")
    all(isapprox.(sum(matrix; dims=2), 1.0; atol=1e-10)) ||
        error("$name transition rows must sum to one")
end
write_table("Table_A.2_DeBacker.tex", debacker_decile_table(DeBackerDataFull))
DeBackerDifference = DeBackerModelFull .- DeBackerDataFull
write_table("Table_A.3_DeBackerDiff.tex",
    debacker_decile_table(DeBackerDifference; differences=true))
# Compute from unrounded probabilities, equally weighting all 100 entries.
DeBackerMAE = 100 * sum(abs, DeBackerDifference) / length(DeBackerDifference)
@printf("%s: mean absolute difference = %.4f percentage points (%.2f rounded).\n",
    "DeBacker comparison",
    DeBackerMAE, DeBackerMAE)
# Relative gaps use each empirical cell probability as the denominator.
# A relative error is undefined if the empirical probability is zero.
if all(x -> x > 0, DeBackerDataFull)
    DeBackerMAPE = 100 * sum(abs, DeBackerDifference ./ DeBackerDataFull) / length(DeBackerDifference)
    @printf("%s: mean absolute relative difference = %.4f%% (%.2f%% rounded).\n",
        "DeBacker comparison",
        DeBackerMAPE, DeBackerMAPE)
else
    println("DeBacker comparison",
        ": mean absolute relative difference is undefined (zero empirical probability).")
end

# Impact decomposition using worker- and entrepreneur-owned investment series.

delta_output = MIT.Output_total_t[1] - SSall.Output_total
abs(delta_output) > eps(Float64) || error("Impact output response is zero")
share_of_output_response(x) = fmt(100x / delta_output, 1)
delta_other = delta_output - (MIT.C_t[1] - SSall.C) - (MIT.Inv_t[1] - SSall.Inv)

io = IOBuffer()
println(io, "\\begin{tabular}{lccccc}")
println(io, "\\toprule")
println(io, " & Cons. & Inv. (total) & Inv. (priv.) & Inv. (market) & Other \\\\")
println(io, "\\midrule")
println(io, "Total & ", share_of_output_response(MIT.C_t[1] - SSall.C), " & ",
    share_of_output_response(MIT.Inv_t[1] - SSall.Inv), " & ",
    share_of_output_response(MIT.InvE_nonM_t[1] - SSall.InvE_nonM), " & ",
    share_of_output_response(MIT.Inv_M_t[1] - SSall.Inv_M), " & ",
    share_of_output_response(delta_other), " \\\\")
println(io, "Workers & ", share_of_output_response(MIT.ConsumpW_t[1] - SSall.ConsumpW), " & ",
    share_of_output_response(MIT.InvCW_t[1] - SSall.InvCW), " & -- & ",
    share_of_output_response(MIT.InvCW_t[1] - SSall.InvCW), " & -- \\\\")
println(io, "Entrepreneurs & ", share_of_output_response(MIT.ConsumpE_t[1] - SSall.ConsumpE), " & ",
    share_of_output_response((MIT.InvE_nonM_t[1] - SSall.InvE_nonM) +
        (MIT.InvCE_t[1] - SSall.InvCE)), " & ",
    share_of_output_response(MIT.InvE_nonM_t[1] - SSall.InvE_nonM), " & ",
    share_of_output_response(MIT.InvCE_t[1] - SSall.InvCE), " & -- \\\\")
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
if CURRENT_OPT.leverage_mode == :riskfree
    write_table("Table_C.3_ShockDecomp.tex", String(take!(io)))
else
    write_table("Table_7_ShockDecomp.tex", String(take!(io)))
end

# Bequest-regime table
bequest_path = joinpath(MODEL_PATH, "bequestregimes.jld2")
if isfile(bequest_path)
    MITbeq_vec, beqsize_vec, beqfrac_vec, beq_stats1 = load(bequest_path,
        "MITbeq_vec", "beqsize_vec", "beqfrac_vec", "beq_stats1")
    all(mit -> mit.converged[1], MITbeq_vec) ||
        error("bequestregimes.jld2 contains an unconverged transition")

    desired_regimes = [(500.0, 1.0), (2000.0, 0.25), (10000.0, 0.05)]
    regimes = map(desired_regimes) do (size, fraction)
        index = findfirst(i -> isapprox(beqsize_vec[i], size) &&
            isapprox(beqfrac_vec[i], fraction), eachindex(beqsize_vec))
        index === nothing && error("Missing bequest regime: size=$size, fraction=$fraction")
        index
    end

    io = IOBuffer()
    println(io, "\\begin{tabular}{lcccccc}")
    println(io, "\\toprule")
    println(io, "& Share with & Bequest size & Top 10\\% & Capital at & Output & Investment \\\\")
    println(io, "Regime & bequest (\\%) & (\\\$1,000) & wealth share & private firms & response & response \\\\")
    println(io, "\\midrule")
    println(io, "Baseline & 0 & 0 & ", fmt(beq_stats1.Top10Share[1], 1), "\\% & ",
        fmt(100beq_stats1.PrivateFirmCapShare[1], 0), "\\% & ",
        fmt(beq_stats1.Y[1], 2), "\\% & ", fmt(beq_stats1.Inv[1], 2), "\\% \\\\")
    println(io, "\\midrule")
    for (number, regime) in enumerate(regimes)
        row = regime + 1
        println(io, number, " & ", fmt(100beq_stats1.ShareBeq[row], 0), " & ",
            replace(fmt(beq_stats1.BeqSize[row], 0), r"000$" => ",000"), " & ",
            @sprintf("%+.1f", beq_stats1.Top10Share[row] - beq_stats1.Top10Share[1]), "pp. & ",
            @sprintf("%+.0f", 100(beq_stats1.PrivateFirmCapShare[row] - beq_stats1.PrivateFirmCapShare[1])), "pp. & ",
            @sprintf("%+.2f", 100(beq_stats1.Y[row] / beq_stats1.Y[1] - 1)), "\\% & ",
            @sprintf("%+.2f", 100(beq_stats1.Inv[row] / beq_stats1.Inv[1] - 1)), "\\% \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write_table("Table_8_BequestRegimes.tex", String(take!(io)))
else
    @warn "Skipping bequest table; result file does not exist" bequest_path
end

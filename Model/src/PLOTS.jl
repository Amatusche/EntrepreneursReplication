"Create the model figures used by the current paper draft."

using CSV, ColorSchemes, Colors, DataFrames, JLD2, LaTeXStrings, PGFPlotsX, StatsPlots
using Plots

include("HANKModel.jl")
using .HANKModel
include("plotPC.jl")
include("plotoutputdec.jl")
pgfplotsx()

const CURRENT_OPT = options()
const MODEL_PATH = CURRENT_OPT.outpath
const DATA_PATH = "../../Data/out/"
const FIGURE_PATH = MODEL_PATH * "Figures/"
mkpath(FIGURE_PATH)

SSall, SSStats, P, G = load(joinpath(MODEL_PATH, "SS.jld2"),
    "SSall", "SSStats", "P", "G")

MIT = load(joinpath(MODEL_PATH, "MIT_Base.jld2"), "MIT")
MIT_entre_exp3 = load(joinpath(MODEL_PATH, "MIT_entre_exp3.jld2"), "MIT")
MIT_worker_exp3 = load(joinpath(MODEL_PATH, "MIT_worker_exp3.jld2"), "MIT")
MITStats = load(joinpath(MODEL_PATH, "MITStats_Base.jld2"), "MITStats")

backcolor = RGB(1, 1, 1)
mypalette = ColorScheme([
    RGB(0.0, 100 / 255, 142 / 255),
    RGB(0.0, 149 / 255, 97 / 255),
    RGB(227 / 255, 41 / 255, 38 / 255),
    RGB(255 / 255, 218 / 255, 71 / 255),
])
theme(:default)
default(titlefont=12, legendfontsize=8, guidefont=10, yminorgrid=false,
    background_color=backcolor, foreground_color_legend=nothing,
    background_color_legend=nothing, linewidth=2, legend=:topright,
    size=(400, 300), legendfonthalign=:left,
    palette=mypalette)

quarters = 8
dt_cum = G.dt_cum .- G.dt[1]
H = findfirst(>=(quarters), dt_cum)

function rates_plot(MIT, SSall, G, dt_cum, H, quarters)
    p = plot(dt_cum[1:H], G.mp_shock_t[1:H] .* 400,
        label=L"Shock ($\epsilon_t$)", xlabel="Quarter",
        ylabel="Percentage points", xlims=(0, quarters), color=:black,
        linewidth=0.5)
    plot!(p, dt_cum[1:H], MIT.infl_t[1:H] .* 400,
        label=L"Inflation ($\pi_t$)", color=1)
    plot!(p, dt_cum[1:H], (MIT.rb_t[1:H] .- SSall.rb) .* 400,
        label=L"Liq. rate ($r^b_t$)", color=1, linestyle=:dot)
    plot!(p, dt_cum[1:H], (MIT.ra_t[1:H] .- SSall.ra) .* 400,
        label=L"Illiq. rate ($r^a_t$)", color=1, linestyle=:dash)
    ylims!(p, -1, 3.0)
    return p
end

# Aggregate output and rate responses: baseline and the two non-response experiments.
poutput = plotoutputdec(MIT, SSall, quarters, H, dt_cum)
poutput_e3 = plotoutputdec(MIT_entre_exp3, SSall, quarters, H, dt_cum)
poutput_w3 = plotoutputdec(MIT_worker_exp3, SSall, quarters, H, dt_cum)
plot!(poutput, title="Baseline", legend=false)
plot!(poutput_e3, title=L"Entrepreneurs do not \\respond to $r^b$", legend=:topright)
plot!(poutput_w3, title=L"Workers do not \\respond to $r^b$", legend=false)
poutput_all = plot(poutput, poutput_e3, poutput_w3, layout=(1, 3), size=(600, 200))
if CURRENT_OPT.leverage_mode == :riskfree
    savefig(poutput_all, FIGURE_PATH * "Figure_C.2_OutputDec_all3.tex")
else
    savefig(poutput_all, FIGURE_PATH * "Figure_4_OutputDec_all3.tex")
end

prates = rates_plot(MIT, SSall, G, dt_cum, H, quarters)
prates_e3 = rates_plot(MIT_entre_exp3, SSall, G, dt_cum, H, quarters)
prates_w3 = rates_plot(MIT_worker_exp3, SSall, G, dt_cum, H, quarters)
plot!(prates, title="Baseline", legend=false)
plot!(prates_e3, title=L"Entrepreneurs do not \\respond to $r^b$", legend=:topright)
plot!(prates_w3, title=L"Workers do not \\respond to $r^b$", legend=false)
prates_all = plot(prates, prates_e3, prates_w3, layout=(1, 3), size=(600, 200))
savefig(prates_all, FIGURE_PATH * "Figure_5_ShockRates_all3.tex")

# Median private-business returns by net-worth decile.
medianreturns = CSV.read(joinpath(DATA_PATH, "medianreturnNWdec.txt"), DataFrame)
return_data = permutedims(Matrix(medianreturns[:, 2:2:(end - 2)]))
preturns = plot(1:10, return_data[1:10, 2], markershape=:cross,
    markerstrokewidth=2, mscolor=1, markeralpha=1, label="SCF 2019")
plot!(preturns, 1:10, return_data[1:10, 1], markershape=:rect,
    markerstrokewidth=2, mscolor=2, markeralpha=1, label="All waves (SCF)")
plot!(preturns, 1:10, SSStats.MedReturn[:, 3], markershape=:circle,
    markerstrokewidth=2, mscolor=3, markeralpha=1, label="Model")
plot!(preturns, xlims=(1, 10), ylims=(0, 1.9), xlabel="Wealth decile",
    ylabel="Median business return", size=(400, 300))
savefig(preturns, FIGURE_PATH * "Figure_3_MedianReturnNW.tex")

# Investment responses and capital ownership in common net-worth groups.
# Markers are ratios of group investment changes to group capital, rather than
# unweighted averages of household ratios. The upper group includes the full tail.
function grouped_investment_ownership(model_path, SSall, P, G)
    Pol, PolE = load(joinpath(model_path, "SS.jld2"), "Pol", "PolE")
    TransPolW, TransPolE = load(joinpath(model_path, "MIT_policies_Base.jld2"),
        "TransPolW", "TransPolE")
    Partial = load(joinpath(model_path, "MITStats_Base.jld2"), "Partial")
    conversion = P.NomGDPFac / 1e6
    nw_e = vec(G.b_grid_MAT_e + G.a_grid_MAT_e + G.e_grid_MAT) .* conversion
    nw_w = vec(G.b_grid_MAT + G.a_grid_MAT) .* conversion
    capital_e = vec(PolE.θ) .* vec(G.e_grid_MAT)
    capital_w = vec(G.a_grid_MAT)
    mass_e = vec(SSall.ge) .* P.share_entre
    mass_w = vec(SSall.g) .* (1 - P.share_entre)
    direct_e = vec(Partial.f_ela_MAT) .* vec(G.e_grid_MAT_adj)
    direct_w = vec(Partial.d_ela_MAT) .* vec(G.a_grid_MAT_adj)
    total_e = vec(TransPolE[1].f) .- vec(PolE.f)
    total_w = vec(TransPolW[1].d) .- vec(Pol.d)
    bounds = [-Inf, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, Inf]
    labels = ["<1", "1-2", "2-4", "4-6", "6-8", "8-10", "10-12", "12-14", "14-16", ">16"]
    function summarize(nw, capital, mass, direct, total)
        capital_mass = capital .* mass
        group_capital = [sum(capital_mass[(bounds[i] .<= nw) .& (nw .< bounds[i+1])])
            for i in 1:length(bounds)-1]
        direct_change = [sum((direct .* mass)[(bounds[i] .<= nw) .& (nw .< bounds[i+1])])
            for i in 1:length(bounds)-1]
        total_change = [sum((total .* mass)[(bounds[i] .<= nw) .& (nw .< bounds[i+1])])
            for i in 1:length(bounds)-1]
        shares = 100 .* group_capital ./ sum(capital_mass)
        return (shares=shares, direct=100 .* direct_change ./ group_capital,
            total=100 .* total_change ./ group_capital)
    end
    return (labels=labels, entrepreneur=summarize(nw_e, capital_e, mass_e, direct_e, total_e),
        worker=summarize(nw_w, capital_w, mass_w, direct_w, total_w))
end
ownership_responses = grouped_investment_ownership(MODEL_PATH, SSall, P, G)
response_values = vcat(ownership_responses.entrepreneur.direct, ownership_responses.entrepreneur.total,
    ownership_responses.worker.direct, ownership_responses.worker.total)
response_min = min(0.0, 0.05 * floor(minimum(response_values) / 0.05))
response_max = max(0.05, 0.05 * ceil(maximum(response_values) / 0.05))
share_max = 8 * ceil(1.15 * maximum(vcat(ownership_responses.entrepreneur.shares,
    ownership_responses.worker.shares)) / 10)

function ownership_response_plot(group, labels, panel_title;
    show_legend=false, show_response_label=false, show_share_label=false)
    x = collect(eachindex(labels))
    # Put bars on the right axis and draw the response markers above them
    # on a mirrored twin axis at the left.
    p = bar(x, group.shares, xticks=(x, labels), xrotation=35,
        xlabel="Net worth (million \$)",
        ylabel=show_share_label ? "Capital share (%)" : "",
        ymirror=true, yguideposition=:right,
        ylims=(0, share_max), xlims=(0.4, length(labels) + 0.6),
        color=:lightgrey, linecolor=:white, bar_width=0.65,
        legend=false, title=panel_title, titlefontsize=10, xgrid=false, ygrid=false)
    response_axis = twinx(p)
    plot!(response_axis, x, group.direct, color=1, marker=:xcross,
        markersize=4, markerstrokewidth=1.5, markercolor=1, markerstrokecolor=1, linewidth=1.5,
        label="Direct response", xlims=(0.4, length(labels) + 0.6),
        ymirror=false, yguideposition=:left,
        ylims=(response_min, response_max),
        ylabel=show_response_label ? "Investment response (%)" : "",
        ygrid=false, legend=show_legend ? :topright : false,
        legend_columns=1)
    plot!(response_axis, x, group.total, color=2, marker=:xcross,
        markersize=4, markerstrokewidth=1.5, markercolor=2, markerstrokecolor=2, linewidth=1.5,
        linestyle=:dash, label="Total response")
    return p
end
pentre = ownership_response_plot(ownership_responses.entrepreneur, ownership_responses.labels,
    "Entrepreneurs: private capital"; show_response_label=true)
pworker = ownership_response_plot(ownership_responses.worker, ownership_responses.labels,
    "Workers: corporate capital"; show_legend=true, show_share_label=true)
pelasticities = plot(pentre, pworker, layout=@layout([a{0.49w} _ b{0.49w}]), size=(600, 265))
if CURRENT_OPT.leverage_mode == :riskfree
    savefig(pelasticities, FIGURE_PATH * "Figure_C.3_ElaFA.tex")
else
    savefig(pelasticities, FIGURE_PATH * "Figure_6_ElaFA.tex")
end

# Bequest-regime figures
bequest_path = joinpath(MODEL_PATH, "bequestregimes.jld2")
if isfile(bequest_path)
    MITStatsbeq_vec, beqsize_vec, beqfrac_vec, beq_stats2 =
        load(bequest_path, "MITStatsbeq_vec", "beqsize_vec",
            "beqfrac_vec", "beq_stats2")

    desired_regimes = [(500.0, 1.0), (2000.0, 0.25), (10000.0, 0.05)]
    regimes = [findfirst(i -> beqsize_vec[i] == size && beqfrac_vec[i] == fraction,
        eachindex(beqsize_vec)) for (size, fraction) in desired_regimes]

    baseline_capital_shares = MITStats.PrivateFirmCapByNWShare
    regime_capital_shares = [MITStatsbeq_vec[i].PrivateFirmCapByNWShare for i in regimes]
    capital_share_changes = 100 .* (hcat(regime_capital_shares...) .-
        baseline_capital_shares)

    group_labels = [raw"<\$2m", raw"\$2m-6m", raw"\$6m-10m", raw">\$10m"]
    regime_labels = [raw"Regime 1 (100\%)", raw"Regime 2 (25\%)", raw"Regime 3 (5\%)"]
    regime_colors = [2, 4, 3]
    regime_patterns = ["north east lines", "crosshatch", "dots"]
    pbeqcap = plot(xticks=(1:4, group_labels),
        ylabel="Change in capital share (pp)", xlabel="Net worth",
        legend=:outertop, legend_columns=3, size=(600, 300), xgrid=false,
        xlims=(0.5, 4.5),
        ylims=(-1, 1) .* (0.5 * ceil(1.15 * maximum(abs, capital_share_changes) / 0.5)))
    hline!(pbeqcap, [0], color=:black, linewidth=0.7, label=false)
    for column in eachindex(regimes)
        # Separate shapes apply patterns to every bar and legend entry.
        for group in 1:4
            center = group + (column - 2) * 0.24
            height = capital_share_changes[group, column]
            rectangle = Shape([center - 0.11, center + 0.11,
                center + 0.11, center - 0.11], [0, 0, height, height])
            plot!(pbeqcap, rectangle,
                label=group == 1 ? regime_labels[column] : false,
                fillcolor=regime_colors[column], linecolor=regime_colors[column],
                linewidth=0.5, fillalpha=1.0,
                extra_kwargs=Dict(:series => Dict("postaction" =>
                    "{pattern=$(regime_patterns[column]), pattern color=black!70}")))
        end
    end
    savefig(pbeqcap, FIGURE_PATH * "Figure_9_BeqCapShares.tex")

    regime_names = collect(1:length(regimes))
    alphas = repeat([1.0 0.0 0.7 0.4], length(regimes), 1)
    response_parts = hcat(beq_stats2.Δ∂IPo∂rdroI[regimes],
        (beq_stats2.ΔIoI[regimes] .- beq_stats2.Δ∂Io∂rdroI[regimes]),
        beq_stats2.Δ∂ICEo∂rdroI[regimes],
        beq_stats2.Δ∂ICWo∂rdroI[regimes]) .* 100 ./ beq_stats2.dIoI
    pbeqresponse = groupedbar(regime_names, response_parts,
        label=["direct private" "general equilibrium" "direct corp. (entrep.)" "direct corp. (workers)"],
        orientation=:h, legend=:outertop, linewidth=0.5, legend_columns=2,
        yflip=true, color=1, fillalpha=alphas)
    plot!(pbeqresponse, yticks=regime_names, ytick_direction=:none,
        ylabel="Regime", ygrid=false, ylims=(0.5, 3.5))
    savefig(pbeqresponse, FIGURE_PATH * "Figure_8a_BeqInvestmentResponse.tex")

    direct_parts = hcat(beq_stats2.ΔϵPIPoI[regimes],
        beq_stats2.ϵPΔIPoI[regimes]) .* 100 ./ beq_stats2.dIoI
    pbeqdirect = groupedbar(regime_names, direct_parts,
        label=["average elasticity" "capital share"], orientation=:h,
        bar_position=:stack, color=[1 1], fillalpha=[0.5 1.0],
        legend=:outertop, linewidth=0.5,
        legend_columns=-1, yflip=true)
    bar_series = filter(series -> series[:seriestype] == :shape, pbeqdirect.series_list)
    for (series, pattern) in zip(bar_series, ("north west lines", "crosshatch"))
        series[:extra_kwargs] = Dict("postaction" =>
            "{pattern=$pattern, pattern color=black!70}")
    end
    plot!(pbeqdirect, yticks=regime_names, ytick_direction=:none,
        ylims=(0.5, 3.5))
    savefig(pbeqdirect, FIGURE_PATH * "Figure_8b_BeqInvestmentResponse_direct.tex")
else
    @warn "Skipping bequest figures; result file does not exist" bequest_path
end

# Phillips-curve appendix figure. The underlying shock-grid experiments are optional.
pc_path = joinpath(MODEL_PATH, "PC.jld2")
pc_bequest_path = joinpath(MODEL_PATH, "PC_beq3.jld2")
if isfile(pc_path) && isfile(pc_bequest_path)
    MITVec = load(pc_path, "MITVec")
    MITbeq_PC = load(pc_bequest_path, "MITVec")
    PhillipsCurve, MCInflatCurve, OutputMCCurve, rbInflCurve = plotPC(MITVec, 0, G, "Baseline")
    plotPC!(PhillipsCurve, MCInflatCurve, OutputMCCurve, rbInflCurve,
        MITbeq_PC, 0, G, "Bequests: 5%, USD 10m", 3, :x)
    plot!(MCInflatCurve, legend=false)
    plot!(OutputMCCurve, legend=false)
    plot!(rbInflCurve, legend=false)
    pc_all = plot(PhillipsCurve, MCInflatCurve, OutputMCCurve, rbInflCurve,
        layout=(2, 2), size=(600, 400))
    savefig(pc_all, FIGURE_PATH * "Figure_C.1_PhillipsCurveAll.tex")
else
    @warn "Skipping Phillips-curve figure; shock-grid results are incomplete" pc_path pc_bequest_path
end

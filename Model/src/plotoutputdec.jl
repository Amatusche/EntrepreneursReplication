"Plot the output decomposition used in the paper."
function plotoutputdec(MIT, SSall::HANKModel.SSstruct, quarters, H, dt_cum)
    poutput = plot(dt_cum[1:H],
        MIT.Output_total_t[1:H] ./ SSall.Output_total .* 100 .- 100,
        title="Output", xlabel="Quarter", ylabel="%", legend=false,
        label=false, yguidefontrotation=-90, xlims=(0, quarters),
        color=:black, linewidth=1)

    lower = zeros(H)
    upper = (MIT.ConsumpW_t[1:H] .- SSall.ConsumpW) ./ SSall.Output_total .* 100
    plot!(poutput, dt_cum[1:H], lower, fillrange=upper, color=1,
        fillalpha=0.5, label="Consumption", linewidth=0, alpha=1,
        legend=:topright)

    lower = upper
    upper = lower .+ (MIT.ConsumpE_t[1:H] .- SSall.ConsumpE) ./ SSall.Output_total .* 100
    plot!(poutput, dt_cum[1:H], lower, fillrange=upper, color=1,
        fillalpha=0.9, label=false, linewidth=0, alpha=1)

    lower = upper
    upper = lower .+ (MIT.InvCW_t[1:H] .- SSall.InvCW) ./ SSall.Output_total .* 100
    plot!(poutput, dt_cum[1:H], lower, fillrange=upper, color=2,
        fillalpha=0.5, label="Market investm.", linewidth=0, alpha=1)

    lower = upper
    upper = lower .+ (MIT.InvCE_t[1:H] .- SSall.InvCE) ./ SSall.Output_total .* 100
    plot!(poutput, dt_cum[1:H], lower, fillrange=upper, color=2,
        fillalpha=0.9, label=false, linewidth=0, alpha=1)

    lower = upper
    upper = lower .+ (MIT.InvE_nonM_t[1:H] .- SSall.InvE_nonM) ./ SSall.Output_total .* 100
    plot!(poutput, dt_cum[1:H], lower, fillrange=upper, color=3,
        fillalpha=0.9, label="Priv. investm.", linewidth=0, alpha=1)

    lower = upper
    upper = MIT.Output_total_t[1:H] ./ SSall.Output_total .* 100 .- 100
    plot!(poutput, dt_cum[1:H], lower, fillrange=upper, color=:gray,
        fillalpha=0.5, label="Other", linewidth=0, alpha=1,
        ylims=(-0.01, 1.8))
    return poutput
end

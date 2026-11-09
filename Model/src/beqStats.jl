function beqStats(SSall, SSStats, MIT, MITStats, Partial, SSallbeq_vec, SSStatsbeq_vec, MITbeq_vec, MITStatsbeq_vec, Partialbeq_vec, beqsize_vec, beqfrac_vec)

    beq_stats = Array{Float64}(undef, length(beqsize_vec) + 1, 20)
    beq_stats[2:end, 1] = beqfrac_vec
    beq_stats[2:end, 2] = beqsize_vec
    for ii = eachindex(beqsize_vec)
        beq_stats[ii+1, 3] = (MITbeq_vec[ii].Output_total_t[1] ./ SSallbeq_vec[ii].Output_total .- 1.0) .* 100.0
        beq_stats[ii+1, 4] = (MITbeq_vec[ii].Inv_t[1] ./ SSallbeq_vec[ii].Inv .- 1.0) .* 100.0
        beq_stats[ii+1, 5] = SSStatsbeq_vec[ii].LaborShareE
        beq_stats[ii+1, 6] = SSStatsbeq_vec[ii].PrivateFirmCapShare
        beq_stats[ii+1, 7] = SSStatsbeq_vec[ii].NWShares[3, 2] .* 100.0
        beq_stats[ii+1, 8] = (MITbeq_vec[ii].rb_t[1] .- SSallbeq_vec[ii].rb) .* 100.0
        beq_stats[ii+1, 9] = MITbeq_vec[ii].infl_t[1] .* 100.0
        beq_stats[ii+1, 10] = MITStatsbeq_vec[ii].f_partial_ela_avg .* 100.0
        beq_stats[ii+1, 11] = MITStatsbeq_vec[ii].f_ela_avg .* 100.0
        beq_stats[ii+1, 12] = beq_stats[ii+1, 3] ./ beq_stats[ii+1, 7]
        beq_stats[ii+1, 13] = beq_stats[ii+1, 3] ./ beq_stats[ii+1, 6] ./ 10.0
        beq_stats[ii+1, 14] = MITStatsbeq_vec[ii].f_partial_ela_avg .* 100.0 ./ beq_stats[ii+1, 6]
        beq_stats[ii+1, 15] = MITStatsbeq_vec[ii].f_ela_avg .* 100.0 ./ beq_stats[ii+1, 6] ./ 10.0
        beq_stats[ii+1, 16] = MITStatsbeq_vec[ii].AggrPartialSemiElaInvE_nonM
        beq_stats[ii+1, 17] = MITStatsbeq_vec[ii].AggrPartialSemiElaInv_M
        beq_stats[ii+1, 18] = MITStatsbeq_vec[ii].AggrPartialSemiElaInvCW
        beq_stats[ii+1, 19] = MITStatsbeq_vec[ii].AggrPartialSemiElaInvCE
        beq_stats[ii+1, 20] = MITStatsbeq_vec[ii].AggrPartialSemiElaInv
    end
    beq_stats[1, 1] = 0.0
    beq_stats[1, 2] = 0.0
    beq_stats[1, 3] = (MIT.Output_total_t[1] ./ SSall.Output_total .- 1.0) .* 100.0
    beq_stats[1, 4] = (MIT.Inv_t[1] ./ SSall.Inv .- 1.0) .* 100.0
    beq_stats[1, 5] = SSStats.LaborShareE
    beq_stats[1, 6] = SSStats.PrivateFirmCapShare
    beq_stats[1, 7] = SSStats.NWShares[3, 2] .* 100.0
    beq_stats[1, 8] = (MIT.rb_t[1] .- SSall.rb) .* 100.0
    beq_stats[1, 9] = MIT.infl_t[1] .* 100.0
    beq_stats[1, 10] = MITStats.f_partial_ela_avg .* 100.0
    beq_stats[1, 11] = MITStats.f_ela_avg .* 100.0
    beq_stats[1, 12] = beq_stats[1, 3] ./ beq_stats[1, 7]
    beq_stats[1, 13] = beq_stats[1, 3] ./ beq_stats[1, 6] ./ 10.0
    beq_stats[1, 14] = MITStats.f_partial_ela_avg .* 100.0 ./ beq_stats[1, 6]
    beq_stats[1, 15] = MITStats.f_ela_avg .* 100.0 ./ beq_stats[1, 6] ./ 10.0
    beq_stats[1, 16] = MITStats.AggrPartialSemiElaInvE_nonM
    beq_stats[1, 17] = MITStats.AggrPartialSemiElaInv_M
    beq_stats[1, 18] = MITStats.AggrPartialSemiElaInvCW
    beq_stats[1, 19] = MITStats.AggrPartialSemiElaInvCE
    beq_stats[1, 20] = MITStats.AggrPartialSemiElaInv

    beq_stats1 = DataFrame(beq_stats, ["ShareBeq", "BeqSize", "Y", "Inv", "LaborShareE", "PrivateFirmCapShare", "Top10Share", "rb", "infl", "fela_partial", "fela", "ΔY/Δinfl", "ΔY/Δrb", "fela_partial/rb", "fela/rb", "AggrSemiElaInvE_nonM", "AggrSemiElaInv_M", "AggrSemiElaInvCW", "AggrSemiElaInvCE", "AggrSemiElaInv"])

    dIoI = MIT.Inv_t[1]/SSall.Inv - 1.0
    ∂Io∂rdroI = (Partial.MIT_partial.Inv_t[1]-SSall.Inv)/SSall.Inv
    dICoI = MIT.Inv_M_t[1]/SSall.Inv_M - 1.0
    ∂ICo∂rdroI = Partial.MIT_partial.Inv_M_t[1]/SSall.Inv_M - 1.0
    dICWoI = (MIT.InvCW_t[1]-SSall.InvCW)/SSall.Inv
    ∂ICWo∂rdroI = (Partial.MIT_partial.InvCW_t[1]-SSall.InvCW)/SSall.Inv
    dICEoI = (MIT.InvCE_t[1] - SSall.InvCE)/SSall.Inv
    ∂ICEo∂rdroI = (Partial.MIT_partial.InvCE_t[1] - SSall.InvCE)/SSall.Inv
    dIPoI = (MIT.InvE_nonM_t[1] - SSall.InvE_nonM)/SSall.Inv
    ∂IPo∂rdroI = (Partial.MIT_partial.InvE_nonM_t[1]-SSall.InvE_nonM)/SSall.Inv

    Δrb = MIT.rb_t[1] - SSall.rb

    dIoI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    ∂Io∂rdroI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    dICoI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    ∂ICWo∂rdroI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    dICEoI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    ∂ICEo∂rdroI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    dIPoI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    ∂IPo∂rdroI_beq = Vector{Float64}(undef, length(MITbeq_vec))
    ΔIoI = Vector{Float64}(undef, length(MITbeq_vec))
    Δ∂Io∂rdroI = Vector{Float64}(undef, length(MITbeq_vec))
    Δ∂ICWo∂rdroI = Vector{Float64}(undef, length(MITbeq_vec))
    ΔϵCW = Vector{Float64}(undef, length(MITbeq_vec))
    ΔϵCWICWoI = Vector{Float64}(undef, length(MITbeq_vec))
    ΔICW = Vector{Float64}(undef, length(MITbeq_vec))
    ϵCWΔICWoI = Vector{Float64}(undef, length(MITbeq_vec))
    Δ∂ICEo∂rdroI = Vector{Float64}(undef, length(MITbeq_vec))
    ΔϵCE = Vector{Float64}(undef, length(MITbeq_vec))
    ΔϵCEICEoI = Vector{Float64}(undef, length(MITbeq_vec))
    ΔICE = Vector{Float64}(undef, length(MITbeq_vec))
    ϵCEΔICEoI = Vector{Float64}(undef, length(MITbeq_vec))
    Δ∂IPo∂rdroI = Vector{Float64}(undef, length(MITbeq_vec))
    ΔϵP = Vector{Float64}(undef, length(MITbeq_vec))
    ΔϵPIPoI = Vector{Float64}(undef, length(MITbeq_vec))
    ΔIP = Vector{Float64}(undef, length(MITbeq_vec))
    ϵPΔIPoI = Vector{Float64}(undef, length(MITbeq_vec))

    for ii = 1:length(MITbeq_vec)
        dIoI_beq[ii] = MITbeq_vec[ii].Inv_t[1]/SSallbeq_vec[ii].Inv - 1.0
        ∂Io∂rdroI_beq[ii] = (Partialbeq_vec[ii].MIT_partial.Inv_t[1]-SSallbeq_vec[ii].Inv)/SSallbeq_vec[ii].Inv
        dICoI_beq[ii] = (MITbeq_vec[ii].Inv_M_t[1]-SSallbeq_vec[ii].Inv_M)/SSallbeq_vec[ii].Inv
        ∂ICWo∂rdroI_beq[ii] = (Partialbeq_vec[ii].MIT_partial.InvCW_t[1]-SSallbeq_vec[ii].InvCW)/SSallbeq_vec[ii].Inv
        dICEoI_beq[ii] = (MITbeq_vec[ii].InvCE_t[1] - SSallbeq_vec[ii].InvCE)/SSallbeq_vec[ii].Inv
        ∂ICEo∂rdroI_beq[ii] = (Partialbeq_vec[ii].MIT_partial.InvCE_t[1]-SSallbeq_vec[ii].InvCE)/SSallbeq_vec[ii].Inv
        dIPoI_beq[ii] = (MITbeq_vec[ii].InvE_nonM_t[1] - SSallbeq_vec[ii].InvE_nonM)/SSallbeq_vec[ii].Inv
        ∂IPo∂rdroI_beq[ii] = (Partialbeq_vec[ii].MIT_partial.InvE_nonM_t[1]-SSallbeq_vec[ii].InvE_nonM)/SSallbeq_vec[ii].Inv

        ΔIoI[ii] = dIoI_beq[ii] - dIoI
        Δ∂Io∂rdroI[ii] = ∂Io∂rdroI_beq[ii] - ∂Io∂rdroI

        # Corporate sector: workers
        # Change in partial response
        Δ∂ICWo∂rdroI[ii] = ∂ICWo∂rdroI_beq[ii] - ∂ICWo∂rdroI
        # Change in average elasticity
        ΔϵCW[ii] = MITStatsbeq_vec[ii].AggrPartialSemiElaInvCW*Δrb-MITStats.AggrPartialSemiElaInvCW*Δrb
        # Total change due to change in average elasticity
        ΔϵCWICWoI[ii] = ΔϵCW[ii] * SSall.InvCW/SSall.Inv
        # Change in capital share
        ΔICW[ii] = SSallbeq_vec[ii].InvCW/SSallbeq_vec[ii].Inv - SSall.InvCW/SSall.Inv
        # Total change due to change in capital share
        ϵCWΔICWoI[ii] = ΔICW[ii]*MITStatsbeq_vec[ii].AggrPartialSemiElaInvCW*Δrb

        # Corporate sector: entrepreneurs
        # Change in partial response
        Δ∂ICEo∂rdroI[ii] = ∂ICEo∂rdroI_beq[ii] - ∂ICEo∂rdroI
        # Change in average elasticity
        ΔϵCE[ii] = MITStatsbeq_vec[ii].AggrPartialSemiElaInvCE*Δrb-MITStats.AggrPartialSemiElaInvCE*Δrb
        # Total change due to change in average elasticity
        ΔϵCEICEoI[ii] = ΔϵCE[ii] * SSall.InvCE/SSall.Inv
        # Change in capital share
        ΔICE[ii] = SSallbeq_vec[ii].InvCE/SSallbeq_vec[ii].Inv - SSall.InvCE/SSall.Inv
        # Total change due to change in capital share
        ϵCEΔICEoI[ii] = ΔICE[ii]*MITStatsbeq_vec[ii].AggrPartialSemiElaInvCE*Δrb

        # Private sector
        # Change in partial response
        Δ∂IPo∂rdroI[ii] = ∂IPo∂rdroI_beq[ii] - ∂IPo∂rdroI
        # Change in average elasticity
        ΔϵP[ii] = MITStatsbeq_vec[ii].AggrPartialSemiElaInvE_nonM*Δrb-MITStats.AggrPartialSemiElaInvE_nonM*Δrb
        # Total change due to change in average elasticity
        ΔϵPIPoI[ii] = ΔϵP[ii] * SSall.InvE_nonM/SSall.Inv
        # Change in capital share
        ΔIP[ii] = SSallbeq_vec[ii].InvE_nonM/SSallbeq_vec[ii].Inv - SSall.InvE_nonM/SSall.Inv
        # Total change due to change in capital share
        ϵPΔIPoI[ii] = ΔIP[ii]*MITStatsbeq_vec[ii].AggrPartialSemiElaInvE_nonM*Δrb
    end

    # ΔϵPIPoI + ϵPΔIPoI  - (∂IPo∂rdroI_beq .- ∂IPo∂rdroI)
    # ΔϵCEICEoI + ϵCEΔICEoI  - (∂ICEo∂rdroI_beq .- ∂ICEo∂rdroI)
    # ΔϵCWICWoI + ϵCWΔICWoI  - (∂ICWo∂rdroI_beq .- ∂ICWo∂rdroI)

    # Δ∂ICWo∂rdroI + Δ∂ICEo∂rdroI + Δ∂IPo∂rdroI

    beq_stats2 = (dIoI=dIoI,
        ∂Io∂rdroI=∂Io∂rdroI,
        dICoI=dICoI,
        ∂ICo∂rdroI=∂ICo∂rdroI,
        dICWoI=dICWoI,
        ∂ICWo∂rdroI=∂ICWo∂rdroI,
        dICEoI=dICEoI,
        ∂ICEo∂rdroI=∂ICEo∂rdroI,
        dIPoI=dIPoI,
        ∂IPo∂rdroI=∂IPo∂rdroI,
        ΔIoI=ΔIoI,
        Δ∂Io∂rdroI=Δ∂Io∂rdroI,
        Δ∂ICWo∂rdroI=Δ∂ICWo∂rdroI,
        ΔϵCW=ΔϵCW,
        ΔϵCWICWoI=ΔϵCWICWoI,
        ΔIC=ΔICW,
        ϵCWΔICWoI=ϵCWΔICWoI,
        Δ∂ICEo∂rdroI=Δ∂ICEo∂rdroI,
        ΔϵCE=ΔϵCE,
        ΔϵCEICEoI=ΔϵCEICEoI,
        ΔICE=ΔICE,
        ϵCEΔICEoI=ϵCEΔICEoI,
        Δ∂IPo∂rdroI=Δ∂IPo∂rdroI,
        ΔϵP=ΔϵP,
        ΔϵPIPoI=ΔϵPIPoI,
        ΔIP=ΔIP,
        ϵPΔIPoI=ϵPΔIPoI)

    return beq_stats1, beq_stats2

end

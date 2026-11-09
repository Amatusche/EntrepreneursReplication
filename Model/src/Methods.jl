module Methods

using SpecialFunctions: erf

export rouwenhorst, tauchen, normalcdf, rand_num_gen, draw_discrete_states,
    mc_simul, own_hist, val_fct_iter, val_fct_update, endog_grid_update, endog_grid_iter,
    lin_interp, gini, indiv_simul, compcrit, discr_norm, bisection, bisection_intervalwidth

"""
    rouwenhorst(Variance, Autocorr, Skew, NumberStates)

Rouwenhorst Method for AR(1) with moments `Expect, Variance, Autocorr, Skew`
"""
function rouwenhorst(Expect, Variance, Autocorr, Skew, NumberStates)
    if Expect == 0 && Skew == 0 # closed form exits
        m = 0
        p = (1 + Autocorr) / 2
        q = p
        Δ = sqrt(Variance) * sqrt(NumberStates - 1)
    else
        # Find parameters numerically
        function moments_diff(x::Array{Float64,2})
            p = x[1, 1]
            q = x[1, 2]
            Δ = x[1, 3]
            m = x[1, 4]
            ExpectDiff = Expect - (m + Δ * (q - p) / (2 - p - q))
            AutocorrDiff = Autocorr - (p + q - 1)
            VarDiff = Variance - (4 * Δ^2 / (NumberStates - 1) * (1 - p) * (1 - q) / ((2 - p - q)^2))
            if (1 - p) * (1 - q) < 0
                SkewDiff = 9999
            else
                SkewDiff = Skew - (p - q) / sqrt((NumberStates - 1) * (1 - p) * (1 - q))
            end
            MomentsDiff = (ExpectDiff)^2 + (VarDiff)^2 + (AutocorrDiff)^2 + (SkewDiff)^2
            # ExpectedDiff, AutocorrDiff, VarDiff, SkewDiffSq
        end
        Init = [Autocorr Autocorr 0.5 0.5]
        OptimRes = minimizer(optimize(moments_diff, Init))
        p = OptimRes[1, 1]
        q = OptimRes[1, 2]
        Δ = OptimRes[1, 3]
        m = OptimRes[1, 4]
    end
    States = collect(linspace(m - Δ, m + Δ, NumberStates))

    # Compute Transition Matrix
    TransMatrix = [p 1-p; 1-q q]
    for nn = 3:NumberStates
        TransMatrix_old = copy(TransMatrix)
        TransMatrix = p * [TransMatrix_old zeros(nn - 1, 1); zeros(1, nn)] +
                      (1 - p) * [zeros(nn - 1, 1) TransMatrix_old; zeros(1, nn)] +
                      (1 - q) * [zeros(1, nn); TransMatrix_old zeros(nn - 1, 1)] +
                      q * [zeros(1, nn); zeros(nn - 1, 1) TransMatrix_old]
        for zz = 2:nn-1
            TransMatrix[zz, :] = TransMatrix[zz, :] / 2
        end
    end
    TransMatrix, States
end

"""
    tauchen(σ2ϵ, ρ, NumberStates)

Tauchen's method to approximate AR(1) pt+1 = rho*pt + eps, eps~N(0,sigmasqr_eps) with Markov chain.
"""
function tauchen(σ2ϵ, ρ, NumberStates)
    States = collect(range(-3 * sqrt(σ2ϵ / (1 - ρ^2)), stop=3 * sqrt(σ2ϵ / (1 - ρ^2)), length=NumberStates))
    d = (States[NumberStates] - States[1]) / (NumberStates - 1)
    TransMatrix = zeros(NumberStates, NumberStates)
    for jj = 1:NumberStates
        for kk = 2:NumberStates-1
            TransMatrix[jj, kk] = normalcdf((States[kk] + d / 2 - ρ * States[jj]) / sqrt(σ2ϵ)) - normalcdf((States[kk] - d / 2 - ρ * States[jj]) / sqrt(σ2ϵ))
        end
        TransMatrix[jj, 1] = normalcdf((States[1] + d / 2 - ρ * States[jj]) / sqrt(σ2ϵ))
        TransMatrix[jj, NumberStates] = 1 - normalcdf((States[NumberStates] - d / 2 - ρ * States[jj]) / sqrt(σ2ϵ))
    end
    TransMatrix, States
end

"""
    normalcdf(x)

CDF of normal distribution
"""
function normalcdf(x)
    1 / 2 * (1 + erf(x / sqrt(2)))
end

"""
`discr_norm(μ,σ,N,gridrange)`

Discretize normal distribution.
"""
function discr_norm(μ, σ, N, gridrange)
    x = collect(range(μ .- gridrange .* σ, stop=μ .+ gridrange .* σ, length=N))
    p = zeros(N)
    p[1] = normalcdf(((x[2].+x[1])./2.0 .- μ) ./ σ)
    for ii = 2:N-1
        p[ii] = normalcdf(((x[ii+1].+x[ii])./2.0 .- μ) ./ σ) - normalcdf(((x[ii].+x[ii-1])./2.0 .- μ) ./ σ)
    end
    p[end] = 1 - normalcdf(((x[end].+x[end-1])./2.0 .- μ) ./ σ)
    return x, p
end

"""
`bisection(f, a, b, tol=0.0)`
"""
function bisection(f, aa, bb, tol=1e-18, MaxIter=10000)
    fa = f(aa)
    fb = f(bb)
    if fa * fb > 0
        error("f(a) and f(b) must have different signs")
    end
    for ii = 1:MaxIter
        c = (aa + bb) / 2
        fc = f(c)
        if abs(fc) < tol
            return c
        elseif fa * fc < 0
            bb = c
            fb = fc
        else
            aa = c
            fa = fc
        end

    end
    @warn("Maximum number of iterations reached in bisection method")
    return (aa + bb) / 2
end

"""
Locate a sign change. Unlike `bisection`, convergence is based on
the width of the interval. This is useful for functions that are not continuous, but have a sign change.
"""
function bisection_intervalwidth(f, lower, upper; xtol=1e-10, maxiter=64)
    f_lower = f(lower)
    f_upper = f(upper)
    f_lower * f_upper <= 0.0 ||
        error("f must have different signs at the bounds")

    for _ = 1:maxiter
        midpoint = 0.5 * (lower + upper)
        f_midpoint = f(midpoint)
        if f_midpoint == 0.0 || upper - lower <= xtol * max(1.0, abs(midpoint))
            return midpoint
        elseif f_lower * f_midpoint < 0.0
            upper = midpoint
            f_upper = f_midpoint
        else
            lower = midpoint
            f_lower = f_midpoint
        end
    end

    return 0.5 * (lower + upper)
end

"""
    rand_num_gen(F,inv,NumberDraws)

This function generates random numbers with cdf `F`. If `F` is invertible, hand its
inverse to the function and specify `inv = 1`, otherwise numerical roofinding will
be used to invert `F`.
"""
function rand_num_gen(F, inv, NumberDraws)
    P = rand(NumberDraws)
    x = zeros(Float64, NumberDraws)
    if inv == 1
        x = F.(P)
    elseif inv == 0
        for pp = 1:NumberDraws
            f(x) = F(x) - P[pp]
            x[pp] = fzero(f, 0.2)
        end
    end
    x
end

"""
    draw_discrete_states(Distr,NumberDraws)

This function draws from a discrete distribution where the pdf across states is
    given by column vector `Distr`.
"""
function draw_discrete_states(Distr, NumberDraws)
    P = rand(NumberDraws)
    CDF = cumsum(Distr)'
    x = sum(P .>= CDF, 2) + 1
end

"""
    mc_simul(InitState,States,TransMatrix,NumberPeriods)

Simulates Markov Chain with `States` and `TransMatrix` for `NumberPeriods` periods
"""
function mc_simul(IndexInitState, States, TransMatrix, NumberPeriods)
    NumberStates = length(States)
    if NumberStates == 1
        SimulValue = States * ones(Float64, NumberPeriods + 1)
        SimulIndex = ones(Int64, NumberPeriods + 1)
    else
        alphas = rand(1, NumberPeriods)
        CumTransMatrix = cumsum(TransMatrix, dims=2) # cumulated transition matrix
        SimulIndex = ones(Int64, NumberPeriods + 1)
        SimulIndex[1, 1] = IndexInitState
        SimulValue = zeros(Float64, NumberPeriods + 1)
        SimulValue[1, 1] = States[SimulIndex[1, 1]]
        for tt = 2:NumberPeriods+1
            ss = 1
            while alphas[1, tt-1] > CumTransMatrix[SimulIndex[tt-1], ss] && ss < NumberStates
                SimulIndex[tt] = copy(SimulIndex[tt]) + 1
                ss = copy(ss) + 1
            end
            SimulValue[tt] = States[SimulIndex[tt]]
        end
    end
    return SimulValue, SimulIndex
end

"""
    own_hist(NumberBins,Data)

Computes histogram
"""
function own_hist(Data, NumberBins)
    Counts = zeros(NumberBins, 1)
    Bins = zeros(NumberBins, 1)
    Binbounds = collect(linspace(minimum(Data), maximum(Data), NumberBins + 1))
    Counts[1, 1] = sum(Binbounds[1, 1] .<= Data .<= Binbounds[2, 1])
    Bins[1, 1] = (Binbounds[1, 1] + Binbounds[2, 1]) / 2
    for bb = 2:NumberBins
        Counts[bb, 1] = sum(Binbounds[bb, 1] .< Data .<= Binbounds[bb+1, 1])
        Bins[bb, 1] = (Binbounds[bb, 1] + Binbounds[bb+1, 1]) / 2
    end
    Freq = Counts ./ length(Data) # compute frequencies
    Freq, Bins
end

"""
val_fct_update(ValueIn::Array{Float64,2},RMat::Array{Float64},IncTrans::Array{Float64,2},β::Float64)

Performs one value function iteration step
"""
function val_fct_update(ValueIn::Array{Float64,2}, RMat::Array{Float64}, IncTrans::Array{Float64,2}, β::Float64)
    NumberAssetGridpoints, ~, NumberIncGridpoints = size(RMat)
    ContVal = ones(NumberAssetGridpoints, NumberIncGridpoints)
    ValueOut = ones(NumberAssetGridpoints, NumberIncGridpoints)
    PolIndex = ones(Int64, NumberAssetGridpoints, NumberIncGridpoints)
    for zz = 1:NumberIncGridpoints
        ContVal[:, zz] = ValueIn[:, :] * IncTrans[zz, :]
        for aa = 1:NumberAssetGridpoints
            ValueOut[aa, zz], PolIndex[aa, zz] = findmax(RMat[aa, :, zz] + β * ContVal[:, zz])
        end
    end
    ValueOut, PolIndex
end

"""
val_fct_iter(ValueIn::Array{Float64,2},IncTrans::Array{Float64,2},IncGrid::Array{Float64},AssetGrid::Array{Float64,1},u,R::Float64,β::Float64,Tol::Float64)
Value function iteration
"""
function val_fct_iter(ValueIn::Array{Float64,2}, RMat::Array{Float64}, IncTrans::Array{Float64,2}, β::Float64, Tol::Float64)
    ValFct = ValueIn
    PolIndex = ones(ValueIn)
    Crit = 1.0
    while Crit > Tol
        ValueOld = copy(ValFct)
        ValFct, PolIndex = val_fct_update(ValueOld, RMat, IncTrans, β)
        Crit = maximum(abs.(ValFct - ValueOld))
    end
    ValFct, PolIndex
end

"""
    lin_interp(x1,y,x2)

Piecewise linear interpolation of the mapping x1->y1 onto the grid x2.
Linear Extrapolation at values of x2 outside x1.
"""
function lin_interp(x1, y1, x2)
    NumberInterpPoints = length(x2)
    NumberFctPoints = length(x1)
    v = length(x1)
    phi = zeros(NumberFctPoints, 1)
    y2 = zeros(NumberInterpPoints, 1)
    for xx = 1:NumberInterpPoints
        phi[1] = (x1[2] - x2[xx]) / (x1[2] - x1[1]) * (x2[xx] <= x1[2])
        if NumberFctPoints >= 3
            for jj = 2:NumberFctPoints-1
                phi[jj] = (x2[xx] - x1[jj-1]) / (x1[jj] - x1[jj-1]) * (x1[jj-1] < x2[xx] <= x1[jj]) +
                          (x1[jj+1] - x2[xx]) / (x1[jj+1] - x1[jj]) * (x1[jj] < x2[xx] <= x1[jj+1])
            end
        end
        phi[NumberFctPoints] = (x2[xx] - x1[NumberFctPoints-1]) / (x1[NumberFctPoints] - x1[NumberFctPoints-1]) * (x2[xx] > x1[NumberFctPoints-1])
        y2[xx] = (y1'*phi)[1]
        if x2[xx] < x1[1]
            y2[xx] = y1[1] + (y1[2] - y1[1]) / (x1[2] - x1[1]) * (x2[xx] - x1[1])
        elseif x2[xx] > x1[NumberFctPoints]
            y2[xx] = y1[NumberFctPoints] + (y1[NumberFctPoints] - y1[NumberFctPoints-1]) / (x1[NumberFctPoints] - x1[NumberFctPoints-1]) * (x2[xx] - x1[NumberFctPoints])
        end
    end
    y2
end

"""
endog_grid_update(PolIn::Array{Float64,2},IncTrans::Array{Float64,2},IncGrid::Array{Float64},AssetGrid::Array{Float64},u_c::Function,u_c_inv::Function,R::Float64,β::Float64)

Endogenous gridpoints method update step
"""
function endog_grid_update(PolIn::Array{Float64,2}, IncTrans::Array{Float64,2}, IncGrid::Array{Float64}, AssetGrid::Array{Float64}, u_c::Function, u_c_inv::Function, R::Float64, β::Float64)
    NumberAssetGridpoints = length(AssetGrid)
    NumberIncGridpoints = length(IncGrid)
    PolConsFutAss = ones(NumberAssetGridpoints, NumberIncGridpoints)
    OptimAssets = ones(NumberAssetGridpoints, NumberIncGridpoints)
    PolNew = ones(NumberAssetGridpoints, NumberIncGridpoints)
    for zz = 1:NumberIncGridpoints
        for aa = 1:NumberAssetGridpoints
            PolConsFutAss[aa, zz] = u_c_inv(β * R * IncTrans[zz, :]' * u_c.(PolIn[aa, :]))
            OptimAssets[aa, zz] = (PolConsFutAss[aa, zz] + AssetGrid[aa] - IncGrid[zz]) / R
        end
        PolNew[:, zz] = lin_interp(OptimAssets[:, zz], PolConsFutAss[:, zz], AssetGrid)
        for aa = 1:NumberAssetGridpoints
            if AssetGrid[aa] < OptimAssets[1, zz]
                PolNew[aa, zz] = R * AssetGrid[aa] + IncGrid[zz] - AssetGrid[1]
            end
        end
    end
    PolNew
end

"""
endog_grid_iter(PolIn::Array{Float64,2},IncTrans::Array{Float64,2},IncGrid::Array{Float64},AssetGrid::Array{Float64},u_c::Function,u_c_inv::Function,R::Float64,β::Float64)

Endogenous gridpoints method iteration
"""
function endog_grid_iter(ConsPol::Array{Float64,2}, IncTrans::Array{Float64,2}, IncGrid::Array{Float64}, AssetGrid::Array{Float64}, u_c::Function, u_c_inv::Function, R::Float64, β::Float64, Tol::Float64)
    Crit = 1
    NumberAssetGridpoints = length(AssetGrid)
    NumberProductGridpoints = length(IncGrid)
    while Crit > Tol
        PolIn = copy(ConsPol)
        ConsPol = endog_grid_update(PolIn, IncTrans, IncGrid, AssetGrid, u_c, u_c_inv, R, β)
        Crit = maximum(abs.(ConsPol - PolIn))
    end
    AssetPol = -ConsPol + R .* repeat(AssetGrid, outer=(1, NumberProductGridpoints)) + (repeat(IncGrid, outer=(1, NumberAssetGridpoints)))'
    ConsPol, AssetPol
end

"""
gini(density,values)

Gini index
"""
function gini(AssetDensity, AssetGrid::Array{Float64,1})
    AssetDensity = AssetDensity / (sum(AssetDensity))
    ShareOwned = AssetGrid .* AssetDensity / (sum(AssetGrid .* AssetDensity))
    CumShareOwned = cumsum(ShareOwned)
    area = zeros(length(CumShareOwned))
    for bb = 1:length(AssetDensity)-1
        area[bb] = (CumShareOwned[bb+1] - 0.5 * (ShareOwned[bb])) * AssetDensity[bb, 1]
    end
    gini = 2 * (0.5 - sum(area))
end

"""
indiv_simul(SimulInitProductState,ConsPolFctEGM,REffective, IncGrid,AssetGrid,ProductTrans,NumberSimulPeriods)

simulates individual's assets from EGM consumption policy function
"""
function indiv_simul(SimulInitProductState, ConsPolFctEGM, REffective, IncGrid, AssetGrid, ProductTrans, NumberSimulPeriods)
    SimulCons = zeros(NumberSimulPeriods)
    SimulAssets = zeros(NumberSimulPeriods)
    SimulInc, SimulIncIndex = mc_simul(SimulInitProductState, IncGrid, ProductTrans, NumberSimulPeriods - 1)
    SimulCons[1] = ConsPolFctEGM[1, SimulInitProductState]
    SimulAssets[1] = -SimulCons[1] + SimulInc[1]
    for tt = 2:NumberSimulPeriods
        SimulCons[tt] = lin_interp(AssetGrid, ConsPolFctEGM[:, SimulIncIndex[tt]], SimulAssets[tt-1])[1, 1]
        SimulAssets[tt] = REffective * SimulAssets[tt-1] + SimulInc[tt] - SimulCons[tt]
    end
    SimulAssets, SimulIncIndex
end

"""
`crit = compcrit(x, xold)`

Compute the absolute relative difference of two vectors.
"""

function compcrit(x,xold)
    absdiff = 0.0

    diff = 0.0
    N = length(x)
    for ii = 1:N

        diff = x[ii] / xold[ii] - 1.0
        absdiff = absdiff + abs(diff)

    end

    return  absdiff/N

end

end

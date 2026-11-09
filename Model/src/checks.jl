"`checks(SSall, Pol, PolE, G, P, Opt)`

Checks whether the steady-state solution is internally consistent."
function checks(SSall, Pol, PolE, G, P, Opt)

# Government budget
GovBCDelta = SSall.G + SSall.Trans + SSall.rb .* SSall.BondDem +
             SSall.subs_payment_union - SSall.Tunion -
             SSall.tax_rev_entre - SSall.tax_rev_worker

NetRev_worker_MAT = P.tau .* (SSall.w .* G.z_grid_MAT[:] .* Pol.lab .+
    SSall.AvgPerCapitaRealProfitsToW .* G.z_grid_MAT[:] ./ G.AvgLabProduct) .-
    (SSall.TransW[:] .- (1.0 .- P.tau) .* SSall.AvgPerCapitaRealProfitsToW .*
     G.z_grid_MAT[:] ./ G.AvgLabProduct)
LabDemE = SSall.LabDemE
OutputE = SSall.OutputE
ProfitTransferE = (1.0 .- P.tauE) .* SSall.AvgPerCapitaRealProfitsToE .*
                  OutputE ./ (SSall.AggrOutputE ./ P.share_entre)
NetRev_entre_MAT = P.tauE .* SSall.TaxableBusinessIncomeE .-
                   (SSall.TransE .- ProfitTransferE)
Gov_sp = (1.0 .- P.share_entre) .* dot(SSall.g[:], NetRev_worker_MAT[:]) .+
         P.share_entre .* dot(SSall.ge[:], NetRev_entre_MAT[:]) .-
         SSall.rb .* SSall.BondDem + SSall.Tunion - SSall.subs_payment_union
GDelta = Gov_sp .-SSall.G

# Transfers received by both household types: government transfers plus
# after-tax distributions of aggregate profits.
AggregateTransfersW = (1.0 .- P.share_entre) .* dot(SSall.g[:], SSall.TransW[:])
AggregateTransfersE = P.share_entre .* dot(SSall.ge[:], SSall.TransE[:])
AggregateTransfers = AggregateTransfersW + AggregateTransfersE
ExpectedTransfers = SSall.Trans +
                    (1.0 .- P.tau) .* SSall.RealProfitsToW +
                    (1.0 .- P.tauE) .* SSall.RealProfitsToE
DeltaTransfers = AggregateTransfers - ExpectedTransfers

# Tax revenue from entrepreneurs
TaxRevE = P.share_entre .* P.tauE .* dot(SSall.TaxableBusinessIncomeE[:], SSall.ge[:])
DeltaTaxRevE = SSall.tax_rev_entre .- TaxRevE

# Tax revenue from workers
TaxRevW = P.tau .* (SSall.w .* SSall.AggrLabSup .+ SSall.RealProfitsToW)
DeltaTaxRevW = SSall.tax_rev_worker .- TaxRevW
NetRevWorker = TaxRevW - (1.0 .- P.share_entre) .* SSall.Trans

# Taxable household income
TaxableIncomeW_MAT = SSall.w .* G.z_grid_MAT[:] .* Pol.lab .+
                     SSall.AvgPerCapitaRealProfitsToW .* G.z_grid_MAT[:] ./ G.AvgLabProduct
AggrTaxableIncomeW = (1.0 .- P.share_entre) .* dot(TaxableIncomeW_MAT, SSall.g[:])
AggrTaxableIncomeE = P.share_entre .* dot(SSall.TaxableBusinessIncomeE[:], SSall.ge[:])

# Investment identities used in findsteady.jl
InvMCheck = P.delta .* SSall.AggrCapSup
expected_depreciation_rate = P.delta_e +
                             P.lambda_e / (1.0 + P.eta) * Opt.entreriskjump
OperatingIncomeE = (1.0 .- P.tauE) .* SSall.RealGrossProfitsE .+
                   (-P.delta_e + P.tauE .* expected_depreciation_rate) .*
                   G.e_grid_MAT
OperatingLossE = min.(OperatingIncomeE, 0.0)
InvEnonMCheck = P.delta_e .* SSall.AggrFirmSav .+
                P.share_entre .* dot(PolE.f[:] .+ OperatingLossE[:], SSall.ge[:])
DeltaInvM = SSall.Inv_M - InvMCheck
DeltaInvEnonM = SSall.InvE_nonM - InvEnonMCheck
DeltaInv = SSall.Inv - (InvMCheck + InvEnonMCheck)

# Return dampening applied to the illiquid asset a. The current model no
# longer applies the old return-dampening formula to private business capital e.
ltau0a = SSall.ra .* (G.a_max .* 0.999) .^ (1.0 .- Opt.ltau)
ReturnDampeningE = ltau0a .* dot(G.a_grid_MAT_e .^ Opt.ltau, SSall.ge) .* P.share_entre
ReturnDampeningW = ltau0a .* dot(G.a_grid_MAT .^ Opt.ltau, SSall.g) .* (1.0 .- P.share_entre)
ReturnDampening = ReturnDampeningE + ReturnDampeningW

# Resource constraint in steady state. At death, private business capital that
# is not transferred directly to newborn entrepreneurs is physically converted
# into goods at book value. `G.inh_wealth` is the population-weighted business
# capital endowment of newborn entrepreneurs, so the corresponding death flow
# is not a goods-producing liquidation. Return dampening is a resource loss.
GrossBusinessCapitalDeathFlow = P.zeta .* SSall.AggrFirmSav
InheritedBusinessCapitalDeathFlow = P.zeta .* G.inh_wealth
LiquidatedBusinessCapital = SSall.AggrFirmSav .- G.inh_wealth
LiquidationProceeds = P.zeta .* LiquidatedBusinessCapital
ResourceSources = SSall.Output_total + LiquidationProceeds
GoodsUses = SSall.C + SSall.Inv + SSall.G + SSall.AdjCosts +
            SSall.CapAdjE + SSall.PortfAdj + SSall.BorrCosts
ResourceUses = GoodsUses + ReturnDampening
ResourceDelta = ResourceUses ./ ResourceSources .- 1.0

return (government=(GovBCDelta=GovBCDelta, GDelta=GDelta),
        transfers=(aggregate=AggregateTransfers, expected=ExpectedTransfers,
                   delta=DeltaTransfers),
        taxes=(worker=TaxRevW, entrepreneur=TaxRevE,
               delta_worker=DeltaTaxRevW, delta_entrepreneur=DeltaTaxRevE),
        taxable_income=(worker=AggrTaxableIncomeW, entrepreneur=AggrTaxableIncomeE),
        investment=(market=InvMCheck, private=InvEnonMCheck,
                    operating_loss=P.share_entre .* dot(OperatingLossE[:], SSall.ge[:]),
                    delta_market=DeltaInvM, delta_private=DeltaInvEnonM,
                    delta_total=DeltaInv),
        return_dampening=(worker=ReturnDampeningW, entrepreneur=ReturnDampeningE,
                          total=ReturnDampening),
        resource=(output=SSall.Output_total,
                  gross_business_capital_death_flow=GrossBusinessCapitalDeathFlow,
                  inherited_business_capital_death_flow=InheritedBusinessCapitalDeathFlow,
                  liquidated_business_capital=LiquidatedBusinessCapital,
                  liquidation_proceeds=LiquidationProceeds,
                  sources=ResourceSources, goods_uses=GoodsUses,
                  return_dampening=ReturnDampening,
                  uses=ResourceUses, delta=ResourceDelta))

end

# Check value functions, policy functions, and distributions

# Plot the value functions for workers
plot(G.a_grid,SSall.v[1,:,1,1],label="low b, low z1, low z2",title="Value function",xlabel="a",ylabel="Value")
plot!(G.a_grid,SSall.v[1,:,end,end],label="low b, high z1, high z2")
plot!(G.a_grid,SSall.v[end,:,end,end],label="high b, high z1, high z2")
plot!(G.a_grid,SSall.v[end,:,1,1],label="high b, low z1, low z2")

plot(G.b_grid, SSall.v[:,1,1,1], label="low a, low z1, low z2", title="Value function", xlabel="b", ylabel="Value")
plot!(G.b_grid, SSall.v[:,end,1,1], label="high a, low z1, low z2")
plot!(G.b_grid, SSall.v[:,1,end,end], label="low a, high z1, high z2")
plot!(G.b_grid, SSall.v[:,end,end,end], label="high a, high z1, high z2")

# Plot the value functions for entrepreneurs
plot(G.a_grid,SSall.ve[1,:,1,1],label="low b, low e, low y",title="Value function",xlabel="a",ylabel="Value")
plot!(G.a_grid,SSall.ve[1,:,end,end],label="low b, high e, high y")
plot!(G.a_grid,SSall.ve[end,:,end,end],label="high b, high e, high y")
plot!(G.a_grid,SSall.ve[end,:,1,1],label="high b, low e, low y")

plot(G.b_grid_e, SSall.ve[:,1,1,1], label="low a, low e, low y", title="Value function", xlabel="b", ylabel="Value")
plot!(G.b_grid_e, SSall.ve[:,end,1,1], label="high a, low e, low y")
plot!(G.b_grid_e, SSall.ve[:,1,end,end], label="low a, high e, high y")
plot!(G.b_grid_e, SSall.ve[:,end,end,end], label="high a, high e, high y")

plot(G.e_grid, SSall.ve[1,1,:,1], label="low a, low b, low y", title="Value function", xlabel="e", ylabel="Value")
plot!(G.e_grid, SSall.ve[1,end,:,end], label="low a, high b, high y")
plot!(G.e_grid, SSall.ve[end,end,:,end], label="high a, high b, high y")
plot!(G.e_grid, SSall.ve[end,1,:,1], label="high a, low b, low y")

# Plot the policy functions for workers
# maximum a holdings:
polc = reshape(Pol.c,G.N_b,G.N_a,G.N_j,G.N_k)
pold = reshape(Pol.d,G.N_b,G.N_a,G.N_j,G.N_k)
pollab = reshape(Pol.lab,G.N_b,G.N_a,G.N_j,G.N_k)
adjcosts = P.χ0 .* abs.(pold) .+ P.χ1 .* abs.(pold ./ G.a_grid_MAT_adj) .^ P.χ2 .* G.a_grid_MAT_adj
polbdot = (1.0 .- P.tau) .* SSall.w .* G.z_grid_MAT .* pollab .+ (SSall.rb .+ P.kappa * (G.b_grid .< 0.0) .+ P.zeta * SSall.zeta_corr ).* G.b_grid_MAT .+ SSall.TransW .- polc .- pold .- adjcosts
poladrift = pold .+ SSall.ra.*G.a_grid_MAT

maxb = findlast(polbdot[:,end,end,end] .> 0.0)

plot(G.a_grid,polc[1,:,1,1],label="low b, low z1, low z2",title="c",xlabel="a",ylabel="c")
plot!(G.a_grid,polc[1,:,end,end],label="low b, high z1, high z2")
plot!(G.a_grid,polc[end,:,end,end],label="high b, high z1, high z2")
plot!(G.a_grid,polc[end,:,1,1],label="high b, low z1, low z2")

plot(G.b_grid, polc[:,1,1,1], label="low a, low z1, low z2", title="c", xlabel="b", ylabel="c")
plot!(G.b_grid, polc[:,end,1,1], label="high a, low z1, low z2")
plot!(G.b_grid, polc[:,1,end,end], label="low a, high z1, high z2")
plot!(G.b_grid, polc[:,end,end,end], label="high a, high z1, high z2")

plot(G.b_grid, polbdot[:,1,1,1], label="low a, low z1, low z2", title="b dot", xlabel="b", ylabel="b drift")
plot!(G.b_grid, polbdot[:,1,end,end], label="low a, high z1, high z2")
plot!(G.b_grid, polbdot[:,end,1,1], label="high a, low z1, low z2")
plot!(G.b_grid, polbdot[:,end,end,end], label="high a, high z1, high z2")

plot(G.a_grid, pold[1,:,1,1], label="low b, low z1, low z2", title="d", xlabel="a", ylabel="d")
plot!(G.a_grid, pold[1,:,end,end], label="low b, high z1, high z2")
plot!(G.a_grid, pold[end,:,end,end], label="high b, high z1, high z2")
plot!(G.a_grid, pold[end,:,1,1], label="high b, low z1, low z2")

plot(G.a_grid, poladrift[1,:,1,1], label="low b, low z1, low z2", title="d", xlabel="a", ylabel="a dot")
plot!(G.a_grid, poladrift[1,:,end,end], label="low b, high z1, high z2")
plot!(G.a_grid, poladrift[maxb,:,end,end], label="high b, high z1, high z2")
plot!(G.a_grid, poladrift[maxb,:,1,1], label="high b, low z1, low z2")

# Plot the policy functions for entrepreneurs
polce = reshape(PolE.ce,G.N_be,G.N_a,G.N_e,G.N_y)
plot(G.a_grid,polce[1,:,1,1],label="low b, low e, low y",title="c",xlabel="a",ylabel="c")
plot!(G.a_grid,polce[1,:,end,end],label="low b, high e, high y")
plot!(G.a_grid,polce[end,:,end,end],label="high b, high e, high y")
plot!(G.a_grid,polce[end,:,1,1],label="high b, low e, low y")

plot(G.b_grid_e, polce[:,1,1,1], label="low a, low e, low y", title="c", xlabel="b", ylabel="c")
plot!(G.b_grid_e, polce[:,end,1,1], label="high a, low e, low y")
plot!(G.b_grid_e, polce[:,1,end,end], label="low a, high e, high y")
plot!(G.b_grid_e, polce[:,end,end,end], label="high a, high e, high y")

polde = reshape(PolE.de,G.N_be,G.N_a,G.N_e,G.N_y)
polfe = reshape(PolE.f,G.N_be,G.N_a,G.N_e,G.N_y)
ltau0e = SSall.re[end,end]*(G.e_max*0.999).^(1-Opt.ltau);
re_MAT = permutedims(repeat(SSall.re,1,1,G.N_be,G.N_a),[3,4,1,2]);
edrift = (re_MAT.-P.delta_e).*G.e_grid_MAT .- ltau0e.*(G.e_grid_MAT.^Opt.ltau);
adjcostsd_e = P.χ0 .* abs.(polde) .+ P.χ1 .* abs.(polde ./ G.a_grid_MAT_e_adj) .^ P.χ2 .* G.a_grid_MAT_e_adj
adjcostsf_e = P.χ0 .* abs.(polfe) .+ P.χ1 .* abs.(polfe ./ G.e_grid_MAT_adj) .^ P.χ2 .* G.e_grid_MAT_adj
polbdote = (SSall.rb .+ P.kappa * (G.b_grid_e .< 0.0) .+ P.zeta * SSall.zeta_corr ) .* G.b_grid_MAT_e .+ P.intoliqshare .* edrift .+ SSall.TransE .- polce .- polde .- polfe .- adjcostsd_e - adjcostsf_e

plot(G.b_grid_e, polbdote[:,1,1,1], label="low a, low e, low y", title="b dot", xlabel="b", ylabel="b drift")
plot!(G.b_grid_e, polbdote[:,1,end,end], label="low a, high e, high y")
plot!(G.b_grid_e, polbdote[:,end,1,1], label="high a, low e, low y")
plot!(G.b_grid_e, polbdote[:,end,end,end], label="high a, high e, high y")

plot(G.a_grid, polde[1,:,1,1], label="low b, low e, low y", title="d", xlabel="a", ylabel="d")
plot!(G.a_grid, polde[1,:,end,end], label="low b, high e, high y")
plot!(G.a_grid, polde[end,:,end,end], label="high b, high e, high y")
plot!(G.a_grid, polde[end,:,1,1], label="high b, low e, low y")

plot(G.e_grid, polfe[1,1,:,1], label="low a, low b, low y", title="f", xlabel="e", ylabel="f")
plot!(G.e_grid, polfe[1,end,:,end], label="low a, high b, high y")
plot!(G.e_grid, polfe[end,end,:,end], label="high a, high b, high y")
plot!(G.e_grid, polfe[end,1,:,1], label="high a, low b, low y")

# Plot the asset distribution for workers
plot(G.a_grid, SSall.g_marg_a)
hline!([0.0], label="0.0")

plot(G.b_grid, SSall.g_marg_b)
vline!([G.b_grid[maxb]], label="max b")
hline!([0.0], label="0.0")

# Plot the asset distribution for entrepreneurs
plot(G.a_grid, SSall.ge_marg_a)
vline!([G.a_grid[maxae]], label="max a")
hline!([0.0], label="0.0")
plot(G.b_grid_e, SSall.ge_marg_b)
vline!([G.b_grid_e[maxbe]], label="max b")
hline!([0.0], label="0.0")
plot(G.e_grid, SSall.ge_marg_e)
vline!([G.e_grid[maxee]], label="max e")
hline!([0.0], label="0.0")

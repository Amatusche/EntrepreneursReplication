function MIT_init(SSall,mp_shock_t,Pol,PolE,G,P)
# initial guesses of aggregate state vector sequence
rb_t = SSall.rb .+ mp_shock_t;
q_t = SSall.q.*ones(G.N_t);
mc_t = SSall.pe.*ones(G.N_t);
wagemarkdown_t = zeros(G.N_t);
rk_t = SSall.rk.*ones(G.N_t);
AggrCapSup_t = SSall.AggrCapSup.*ones(G.N_t);
KE_t = SSall.KE.*ones(G.N_t);
KW_t = SSall.KW.*ones(G.N_t);
Inv_t = SSall.Inv.*ones(G.N_t);
Inv_M_t = SSall.Inv_M.*ones(G.N_t);
InvE_t = SSall.InvE.*ones(G.N_t);
InvC_M_t = SSall.InvC_M.*ones(G.N_t);
InvE_M_t = SSall.InvE_M.*ones(G.N_t);
InvE_nonM_t = SSall.InvE_nonM.*ones(G.N_t);
InvCW_t = SSall.InvCW.*ones(G.N_t);
InvCE_t = SSall.InvCE.*ones(G.N_t);

infl_t = P.infl.*ones(G.N_t);
inflw_t = P.inflw.*ones(G.N_t);

AggrLabDemC_t = SSall.AggrLabDemC.*ones(G.N_t);
AggrCapDemC_t = SSall.AggrCapDemC.*ones(G.N_t);
AggrLabSup_t = SSall.AggrLabSup.*ones(G.N_t);
AggrLiqSav_t = SSall.AggrLiqSav.*ones(G.N_t+1);
AggrLiqSavW_t = SSall.AggrLiqSavW.*ones(G.N_t+1);
AggrLiqSavE_t = SSall.AggrLiqSavE.*ones(G.N_t+1);
AggrIlliqSav_t = SSall.AggrIlliqSav.*ones(G.N_t+1);
AggrIlliqSavW_t = SSall.AggrIlliqSavW.*ones(G.N_t+1);
AggrIlliqSavE_t = SSall.AggrIlliqSavE.*ones(G.N_t+1);
G_t = SSall.G.*ones(G.N_t);

ConsumpW_t = SSall.ConsumpW.*ones(G.N_t);
ConsumpE_t = SSall.ConsumpE.*ones(G.N_t);
PortfAdjW_t = SSall.PortfAdjW.*ones(G.N_t);
PortfAdjE_t = SSall.PortfAdjE.*ones(G.N_t);
CapAdjE_t = SSall.CapAdjE.*ones(G.N_t);
BorrCosts_t = SSall.BorrCosts.*ones(G.N_t);
PortfAdj_t = SSall.PortfAdj.*ones(G.N_t);

zeta_corr_t =ones(G.N_t).*SSall.zeta_corr;

if P.share_entre > 0
    AggrOutputE_t = SSall.AggrOutputE.*ones(G.N_t);
    AggrFirmSav_t = SSall.AggrFirmSav.*ones(G.N_t+1);
    AggrLabDemE_t = SSall.AggrLabDemE.*ones(G.N_t);
    AggrCapDemE_t = SSall.AggrCapDemE.*ones(G.N_t);
    tax_rev_entre_t = SSall.tax_rev_entre.*ones(G.N_t);
elseif P.share_entre == 0
    AggrOutputE_t = zeros(G.N_t);
    AggrFirmSav_t = zeros(G.N_t+1);
    AggrLabDemE_t = zeros(G.N_t);
    AggrCapDemE_t = zeros(G.N_t);
    tax_rev_entre_t = zeros(G.N_t);
end
subs_payment_union_t = SSall.subs_payment_union.*ones(G.N_t);

OutputC_t = SSall.OutputC.*ones(G.N_t);
Output_total_t = SSall.Output_total.*ones(G.N_t);
omega_entre = SSall.omega_entre.*ones(G.N_t);
dinfl_t = zeros(G.N_t-1);
dOutput_t = zeros(G.N_t-1);
dinflw_t = zeros(G.N_t-1);
dLab_t = zeros(G.N_t-1);
pe_t = SSall.pe.*ones(G.N_t);
w_t = ones(G.N_t).*SSall.w;
what_t = ones(G.N_t).*SSall.what;
AdjCosts_t = SSall.AdjCosts.*ones(G.N_t);
RealProfits_t = SSall.RealProfits.*ones(G.N_t);
ra_t = SSall.ra.*ones(G.N_t);
Trans_t = SSall.Trans.*ones(G.N_t);
UnionProfits_t = ones(G.N_t).*SSall.UnionProfits;
C_t = SSall.C.*ones(G.N_t);
Excess_Sup_Liq_t = (SSall.AggrLiqSav .- SSall.BondDem).*ones(G.N_t+1);

RealProfitsToW_t= SSall.RealProfitsToW.*ones(G.N_t);
RealProfitsToE_t=  SSall.RealProfitsToE.*ones(G.N_t);
AvgPerCapitaRealProfitsToW_t=  SSall.AvgPerCapitaRealProfitsToW.*ones(G.N_t);
AvgPerCapitaRealProfitsToE_t=  SSall.AvgPerCapitaRealProfitsToE.*ones(G.N_t);
Tunion_t = SSall.Tunion.*ones(G.N_t);

it = [0]

MIT = MITstruct(mp_shock_t, rb_t, q_t, mc_t, wagemarkdown_t, rk_t, AggrCapSup_t,
KW_t, KE_t, Inv_t, Inv_M_t, InvE_t, InvE_nonM_t, InvE_M_t, InvC_M_t, InvCW_t, InvCE_t,
infl_t, inflw_t, AggrCapDemC_t, AggrCapDemE_t, AggrLabDemC_t, AggrLabDemE_t, AggrLabSup_t, AggrLiqSav_t,
AggrLiqSavW_t, AggrLiqSavE_t, AggrIlliqSav_t, AggrIlliqSavW_t, AggrIlliqSavE_t,
ConsumpW_t, ConsumpE_t, PortfAdjW_t, PortfAdjE_t, CapAdjE_t, BorrCosts_t, PortfAdj_t,
AggrOutputE_t, AggrFirmSav_t, tax_rev_entre_t, subs_payment_union_t, Tunion_t, OutputC_t,
Output_total_t, omega_entre, zeta_corr_t, dinfl_t, dOutput_t, dinflw_t, dLab_t, pe_t, w_t, what_t, AdjCosts_t,
RealProfits_t, ra_t, Trans_t, G_t, C_t, Excess_Sup_Liq_t, RealProfitsToW_t,RealProfitsToE_t,
AvgPerCapitaRealProfitsToW_t,AvgPerCapitaRealProfitsToE_t, UnionProfits_t,it,[false],[Inf, Inf])

TransPolW = [deepcopy(Pol) for tt = 1:G.N_t]
TransPolE = [deepcopy(PolE) for tt = 1:G.N_t]

return MIT::MITstruct, TransPolW, TransPolE

end

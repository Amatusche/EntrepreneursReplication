function steady_init(G,P,Opt)
## Construct initial guess for steady state equilibrium

    # Initial guess consists of: v, ve, rb, ra,
    # Workers
    v = zeros(G.N_b,G.N_a,G.N_j,G.N_k);
    lab0 = zeros(G.N_b,G.N_a,G.N_j,G.N_k); # labor supply at zero savings
    # Entrepreneurs
    ve = zeros(G.N_be,G.N_a,G.N_e,G.N_y);

    zeta_corr=1
    rb = 0.005;

    Rb = rb .+ P.kappa*(G.b_grid.<0) .+ P.zeta*zeta_corr;
    AggrLabDemC = 0.445;
    AggrCapSup = 13.24;
    CapLabRatio = AggrCapSup/AggrLabDemC;

    mc = 1.0.-1.0./P.epsilon + P.theta*P.ρ*P.infl/P.epsilon;
    pb = mc;
    pe = mc;
    pp = mc;
    CapOutRatio = CapLabRatio^(1-P.alphaa);
    w = (1-P.alphaa)*mc*(CapLabRatio)^P.alphaa;
    what = (rb*P.inflw)*P.thetaw/P.epsilonw + (1.0 +P.tauW)*(1.0-1.0/P.epsilonw)*w  # wage paid to workers
    Output_public = AggrLabDemC^((1-P.alphaa))*AggrCapSup^(P.alphaa);

    Output_total = 1.699;
    AggrOutputE = 0.3494;
    AdjCosts = P.theta/2*P.infl^2*Output_total;
    rk = P.alphaa*mc*(CapLabRatio)^(P.alphaa-1);
    RealProfits = (1-pb)*Output_total - AdjCosts;
    ra = rk - P.delta +P.zeta;
    Trans = 0.06*Output_total
    TransW = Trans .+ (1/(1-P.share_entre))*(1-P.tau)*(1-P.omega_liq).*G.z_grid_MAT./G.AvgLabProduct.*RealProfits ;
    omega_entre = 0.0

    RealProfitsToW = (1.0 .- P.omega_liq) .* RealProfits # total profits that go to workers liquid account in proportion to labor income
    RealProfitsToE = P.omega_liq .* omega_entre .* RealProfits # total profits that go to entrepreneurs in proportion to output income
    AvgPerCapitaRealProfitsToE = RealProfitsToE ./ P.share_entre # average profits PER ENTREPRENEURS going to business income

    if P.share_entre > 0
        LabDemENoLev = (pe.*P.entre_prod.*reshape(G.y_grid,1,1,1,G.N_y).*(1.0.-P.alphaaE).*P.nu.*G.e_grid_MAT.^(P.alphaaE*P.nu)/w).^(1/(1-(1-P.alphaaE)*P.nu)); # optimal labor demand of entrepreneurs
        OutputE = P.entre_prod*reshape(G.y_grid,1,1,1,G.N_y).*(G.e_grid_MAT.^P.alphaaE.*LabDemENoLev.^(1-P.alphaaE)).^P.nu;
        omega_entre = AggrOutputE/Output_total;
        OutputEnoLev = P.entre_prod .* reshape(G.y_grid,1,1,1,G.N_y) .* (G.e_grid_MAT .^ P.alphaaE .* LabDemENoLev .^ (1.0 .- P.alphaaE)) .^ P.nu
        RealGrossProfitsEnoLev = pe .* OutputEnoLev .- w .* LabDemENoLev
        TransE = Trans .+ (1.0 .- P.tauE) .* AvgPerCapitaRealProfitsToE .* OutputE ./ (AggrOutputE ./ P.share_entre)
        CapDemE = 0.0
    end

    # Initial guesses of value functions
    # Workers
    x0=1.2;
    for j = 1:G.N_j
        for k = 1:G.N_k
            for bb = 1:G.N_b
                lab_eq(x) = lab_fun(x,P.tau,what,G.z_grid[j,k],Rb[bb],G.b_grid[bb],TransW[bb,1,j,k],P.varphi,P.sigma,P.gamma,Opt.min_c);
                lab0[bb,:,j,k] .= min(fzero(lab_eq,x0),P.max_l);
                for aa = 1:G.N_a
                    if P.sigma > 1
                        v[bb,aa,j,k] = ((max.((1.0.-P.tau).*what.*lab0[bb,aa,j,k].*G.z_grid[j,k].+Rb[1].*G.b_grid[bb].+ra.*G.a_grid[aa].+TransW[bb,aa,j,k],Opt.min_c).^(1.0.-P.sigma))./(1.0.-P.sigma))./P.ρ .- (P.varphi.*(lab0[bb,aa,j,k].^(1.0.+1.0./P.gamma))./(1.0.+1.0./P.gamma)) ./ P.ρ;    # initial guess is "to stay put"
                    elseif P.sigma == 1
                        v[bb,aa,j,k] = log.(max.((1.0.-P.tau).*what.*lab0[bb,aa,j,k].*G.z_grid[j,k].+Rb[1].*G.b_grid[bb].+ra.*G.a_grid[aa].+TransW[bb,aa,j,k],Opt.min_c))./P.ρ .- (P.varphi.*(lab0[bb,aa,j,k].^(1.0.+1.0./P.gamma))./(1.0.+1.0./P.gamma)) ./ P.ρ;    # initial guess is "to stay put" ;
                    end
                end
            end
        end
    end

    # Entrepreneurs
    if P.share_entre > 0
    for yy = 1:G.N_y
        for bb = 1:G.N_be
            for aa = 1:G.N_a
                for ee = 1:G.N_e
                    if P.sigma > 1
                        ve[bb,aa,ee,yy] = ((max.(RealGrossProfitsEnoLev[bb,aa,ee,yy].+TransE[bb,aa,ee,yy].+Rb[1].*G.b_grid_e[bb].+ra.*G.a_grid[aa].+TransE[bb,aa,ee,yy],Opt.min_c).^(1.0.-P.sigma))./(1.0.-P.sigma))./P.ρE .- (P.varphi.*(P.entre_lab.^(1.0.+1.0./P.gamma))./(1.0.+1.0./P.gamma)) ./ P.ρE;    # initial guess is "to stay put"
                    elseif P.sigma == 1
                        ve[bb,aa,ee,yy] = log.(max(RealGrossProfitsEnoLev[bb,aa,ee,yy].+TransE[bb,aa,ee,yy].+Rb[1].*G.b_grid_e[bb].+ra.*G.a_grid[aa].+TransE[bb,aa,ee,yy],Opt.min_c))./P.ρE .- (P.varphi.*(P.entre_lab.^(1.0.+1.0./P.gamma))./(1.0.+1.0./P.gamma)) ./ P.ρE;    # initial guess is "to stay put" ;
                    end # use unadjusted entrepreneurial profits for initial guess here, so as to get curvature of ve in e
                end
            end
        end
    end
    end

    gtemp = zeros(G.N);
    gtemp[1] = 1.0;
    gtempE = zeros(G.Ne);
    gtempE[1] = 1.0;

    InitGuess = SSminstruct(v, ve, lab0, TransW, TransE, gtemp, gtempE, rb, ra, pb, pe, pp, w, what, Output_total, AggrOutputE, RealProfits,
    AggrCapSup, CapLabRatio, omega_entre, zeta_corr);

## Find steady state
return InitGuess

end

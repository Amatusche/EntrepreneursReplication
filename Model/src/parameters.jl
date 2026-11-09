function parameters(Opt)

    epsilonw = Inf # elasticity of substitution between labor varieties
    LabShare = 0.67
    alphaa = 1 - LabShare # capital share in corporate sector production
    nu = 0.8 # Published paper of Tan has 0.786
    eta = 2.0 # rate parameter of the exponential distribution of capital destruction shocks (in case of jump risk)
    lambda_e = 1/20

    if Opt.leverage_mode == :none
        calib = [0.016661674759106154, 1.2509583487211346, 0.015648624907548033, 0.3926270706224035, 1.326, 2.075715489928884, 0.9553702051142433, 0.12303690767230138, 0.12870847094814336, 2.0, 0.33, 0.01798, 0.024396, 7.4014, 0.10055258204424791, 0.11668634593888666]
    elseif Opt.leverage_mode == :riskfree
        calib = [0.016669578387353366, 1.0337458907379964, 0.01994501404943988, 0.3995616103011289, 1.326, 2.0081750612817086, 0.5730010452920028, 0.1096226263229056, 0.10231519967191886, 2.0, 0.33, 0.01798, 0.024396, 7.4014, 0.12966318309458166, 0.10350706531772258]
    elseif Opt.leverage_mode == :risky
        calib = [0.016661674759106154, 1.2509583487211346, 0.015648624907548033, 0.3926270706224035, 1.326, 2.075715489928884, 0.9553702051142433, 0.12303690767230138, 0.12870847094814336, 2.0, 0.33, 0.01798, 0.024396, 7.4014, 0.10055258204424791, 0.11668634593888666]
    end
    # ρ               ρE/ρ        κ                    χ1              χ2  entre_prod      log_y_spread   σe             χe1       χe2  αe δe  λe  η  avgtransrate lowstateshare

    sigma_e= calib[8]

    alphaaE = calib[11]

    P = (share_entre=0.071,
        sigma=1, # risk aversion

        # utility function is u(c,l) = c^(1-sigma)/(1-sigma) -varphi*l^(1+1/gamma)/(1+1/gamma)
        zeta=1 / 180, # death rate

        AdjCosts=1,
        ρ=calib[1],
        ρE=calib[1]*calib[2],
        kappa=calib[3],
        χ0=0.0,
        χ1=calib[4],
        χ2=calib[5], # the value in Alves, Kaplan, Moll, Violante is 1.326
        entre_prod=calib[6],

        log_y_spread=calib[7], # spread between the two entrepreneurial productivity states on the log grid
        lambda_e = calib[13], # capital quality shock arrival rate
        eta = calib[14], # parameter of the exponential distribution of capital destruction shocks (in case of jump risk)
        sigma_e=sigma_e, # capital quality shock standard deviation
        χe1=calib[9],
        χe2=calib[10],
        # The productivity states are statistical income states, so both their
        # stationary shares and transition speed are internally calibrated.
        lambda_y_lh=NaN,
        lambda_y_hl=NaN,
        avgtransrate=calib[15], # starting value; calibrated to entrepreneurial-income mobility
        lowstateshare=calib[16],

        varphi=2.243356002397100, # labor disutility factor
        gamma=1.0, # Frisch elasticity
        phi=1.25, # Taylor rule

        underline_b= -1.0, # borrowing constraint
        underline_be = -1.0, # borrowing constraint for entrepreneurs
        alphaa=alphaa, # Capital share
        epsilon=10.0, # Governs demand elasticity final goods sector
        theta=100.0,
        epsilonw=epsilonw, # elasticity of substitution between labor types
        thetaw=0.0, # adjustment cost parameter for wages
        delta= 1-(1- 0.07)^(1/4),
        delta_e = calib[12]- (1 - calib[14] / (calib[14] + 1)) * calib[13] * Opt.entreriskjump, # depreciation net of expected capital destruction
        tau=0.3, # tax on labor income, 0.3 in HANK
        tauE=0.3, # tax on entrepreneurial profits
        tauW=0.0,
        max_l=1.0, # upper bound on labor supply
        entre_lab=1 / 2, # own labor supply of entrepreneurs
        nu=nu,
        alphaaE=alphaaE,
        infl=0.0, # price inflation
        inflw=0.0, # wage inflation
        inheritance_index_e=1,
        frac_inh_wealth=0.0,# share of entrepreneurs inheritance
        transshare=0.06, # output share of transfers
        λ = 1.5, # leverage constraint k < λ * equity (Ang, Cole Lawson find median of LR = liabilities/total assets = 0.305 <=> λ = 1/(1-LR) = 1.43)

        # Factor to convert into USD (Output_total is quarterly GDP per capita in the model, US GDP per capita was 59531.66)
        NomGDPFac=(((21372582 * 1e6 / 4) / (128579 * 1e3)) / 1.82), # ((Annual GDP 2019 / 4)/ Total HHs 2019) / SS.Output_total

        # Factor to convert employment
        EmpFac= 158.6*1e6, # number of employed in 2019

        # Share of entrepreneurial income that goes directly into the liquid
        # account
        intoliqshare=1,
        omega_liq=NaN)

    transition_share_term = P.lowstateshare^2 + (1 - P.lowstateshare)^2
    @set! P.lambda_y_lh = P.avgtransrate * transition_share_term / P.lowstateshare
    @set! P.lambda_y_hl = P.avgtransrate * transition_share_term / (1 - P.lowstateshare)

    @set! P.omega_liq = P.alphaa # profit share reinvested into illiquid account

    return P
end

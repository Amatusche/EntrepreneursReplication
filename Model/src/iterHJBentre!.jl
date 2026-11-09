function iterHJBentre!(ve, PolE, Aentre, Bentre, RbE, rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants,
    adriftE, P, G, Opt, TempE::EntreTemp)
    crit = Inf
    for it = 1:Opt.MaxIterHJB
        crit = HHsolve_entre!(ve, PolE, Aentre, Bentre, RbE, rk, TransE, RealGrossProfitsnoLevPerCap, FirmIncomeConstants,
            adriftE, P, G, Opt, Opt.Delta_e, true, true, TempE)
        if Opt.display == 2
            println("     ITERATION = " * string(it) * ", Crit " * string(crit))

        end
        if crit < Opt.Vtol

            return
        end
    end
    if Opt.display > 0
        println("   Entrepeneur HJB did not converge after ", Opt.MaxIterHJB, " iterations, final crit = ", crit)
    end
end

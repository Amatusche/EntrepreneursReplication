function findsteady_fixedB(SSInitGuess, AggrLiqSav, P, G, Opt)
    tol = Opt.tol
    println("Liquid asset demand = "*string(AggrLiqSav))
    SS, SSall, Pol, PolE, A, Aentre = findsteady(SSInitGuess, P, G, Opt)
    println("Liquid asset supply = "*string(SSall.AggrLiqSav))
    diffsav = AggrLiqSav - SSall.AggrLiqSav
    if abs(diffsav) < tol
        println("Steady state found, rb = ",SSall.rb)
        return SS, SSall, Pol, PolE, A, Aentre
    elseif diffsav  < -tol # too much savings, rate has to be lower than SSall.rb
        println("     Excess liquid saving = ", -diffsav, " at rb = ", SSall.rb)
        interval = [0.002, 1.01*SSall.rb]
    elseif diffsav > tol # too little savings, rate has to be higher than SSall.rb
        println("     Excess liquid saving = ", -diffsav, " at rb = ", SSall.rb)
        interval = [0.99*SSall.rb, 0.01]
    end
    g(x) = excess_liq(x,SSall,AggrLiqSav,P,G,Opt)
    if abs(diffsav) < 1e-6
        rb = find_zero(g, SS.rb, atol = tol)
    else
        rb = find_zero(g, interval, atol = tol)
    end
    println("Steady state found, rb = ",rb)
    @set! SSall.rb = rb
    SS, SSallnew, Pol, PolE, A, Aentre = findsteady(SSall, P, G, Opt)
    return SS, SSallnew, Pol, PolE, A, Aentre
end

function excess_liq(rb,SSall,AggrLiqSav,P,G,Opt)
    @set! SSall.rb = rb

    _, SSallnew, _, _, _, _ = findsteady(SSall, P, G, Opt)
    diffsav = AggrLiqSav - SSallnew.AggrLiqSav
    println("     rb = ",SSall.rb,", excess liquid saving = ",-diffsav)
    return diffsav
end

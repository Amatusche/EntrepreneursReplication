function iterHJB!(v, Pol, A, B, Rb, w, Trans, lab0, adrift, P, G, Opt, Temp)
    for it = 1:Opt.MaxIterHJB
        crit = HHsolve!(v, Pol, A, B, Rb, w, Trans, lab0, adrift,
            P, G, Opt, Opt.Delta, Temp)
        if Opt.display == 2
            println("     ITERATION = " * string(it) * ", Crit " * string(crit))

        end
        if crit < Opt.Vtol
            break
        end
    end

end

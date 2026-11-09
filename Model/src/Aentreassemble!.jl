function Aentreassemble!(Aentre,Y,Xa,Xe,Xb,Za,Ze,Zb,G)
    N = length(G.vaexo)

    for kk = 1:N
        Aentre.nzval[kk] = G.vaexo[kk]
        if G.iaexo[kk] == G.jaexo[kk]
            Aentre.nzval[kk] += Y[G.jaexo[kk]]
        end
        if G.iaexo[kk].-G.N_abe == G.jaexo[kk]
            Aentre.nzval[kk] += Xe[G.jaexo[kk]]
        end
        if G.iaexo[kk].-G.N_be == G.jaexo[kk]
            Aentre.nzval[kk] += Xa[G.jaexo[kk]]
        end
        if G.iaexo[kk].-1 == G.jaexo[kk]
            Aentre.nzval[kk] += Xb[G.jaexo[kk]]
        end
        if G.iaexo[kk].+1 == G.jaexo[kk]
            Aentre.nzval[kk] += Zb[G.iaexo[kk]]
        end
        if G.iaexo[kk].+G.N_be == G.jaexo[kk]
            Aentre.nzval[kk] += Za[G.iaexo[kk]]
        end
        if G.iaexo[kk].+G.N_abe == G.jaexo[kk]
            Aentre.nzval[kk] += Ze[G.iaexo[kk]]
        end
    end

end

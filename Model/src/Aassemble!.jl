function Aassemble!(A,Y,Xa,Xb,Za,Zb,G)
    N = length(G.vaz)

    for kk = 1:N
        A.nzval[kk] = G.vaz[kk]
        if G.iaz[kk] == G.jaz[kk]
            A.nzval[kk] += Y[G.jaz[kk]]
        end
        if G.iaz[kk].-G.N_b == G.jaz[kk]
            A.nzval[kk] += Xa[G.jaz[kk]]
        end
        if G.iaz[kk].-1 == G.jaz[kk]
            A.nzval[kk] += Xb[G.jaz[kk]]
        end
        if G.iaz[kk].+1 == G.jaz[kk]
            A.nzval[kk] += Zb[G.iaz[kk]]
        end
        if G.iaz[kk].+G.N_b == G.jaz[kk]
            A.nzval[kk] += Za[G.iaz[kk]]
        end
    end
end

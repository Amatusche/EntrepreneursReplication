function bbuild!(b,Ne,ce,ce_B,ce_F,Ic_F,Ic_B,Ic_0,lab,ve,Delta_e,P)
    u = 0.0
    gammaexp = (1.0 .+ 1.0 ./ P.gamma)

    if P.sigma > 1.0
        for xx = 1:Ne
            ce[xx] = ce_F[xx] .* Ic_F[xx] .+ ce_B[xx] .* Ic_B[xx] .+ ce[xx] .* Ic_0[xx]
            u = (ce[xx] .^ (1.0 - P.sigma)) ./ (1.0 - P.sigma) .-  P.varphi .* lab[xx].^ gammaexp./gammaexp #z[xx]
            b[xx] = u .+ 1.0 ./ Delta_e .* ve[xx]
        end
    elseif P.sigma == 1.0
        for xx = 1:Ne
            if Ic_F[xx]
                ce[xx] = ce_F[xx]
            elseif Ic_B[xx]
                ce[xx] = ce_B[xx]
            end

            u = log.(ce[xx]) .-  P.varphi .* lab[xx].^ gammaexp./gammaexp
            b[xx] = u .+ ve[xx] ./ Delta_e
        end
    end
end

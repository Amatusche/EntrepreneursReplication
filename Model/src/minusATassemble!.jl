function minusATassemble!(minusAT_e,Aentre,x,i,j)

    transpose!(minusAT_e,Aentre)
    for kk in eachindex(minusAT_e.nzval)
        minusAT_e.nzval[kk] = - minusAT_e.nzval[kk]
    end
    for col = 1:size(minusAT_e, 2)
        for kk = minusAT_e.colptr[col]:(minusAT_e.colptr[col+1]-1)
            if minusAT_e.rowval[kk] == col
                minusAT_e.nzval[kk] += x
                break
            end
        end
    end
end

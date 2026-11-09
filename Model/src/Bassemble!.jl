function Bassemble!(B,x,Aentre,iaexo,jaexo)
N = length(iaexo)
for kk = 1:N
    B.nzval[kk] = - Aentre.nzval[kk]
    if iaexo[kk] == jaexo[kk]
        B.nzval[kk] += x
    end
end

end

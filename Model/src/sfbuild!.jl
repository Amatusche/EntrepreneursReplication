function sfbuild!(sf_F,a,b,c,P)
N = length(sf_F)
for xx = 1:N
    sf_F[xx] = -a[xx]
    if a[xx] < 0.0
        sf_F[xx] -= P.c_n ./ 2.0 .* a[xx] .^ 2.0 ./ b[xx] .+ (1.0 .- P.c_lambda) .* c[xx] .* b[xx]
    elseif a[xx] > 0.0
        sf_F[xx] -= P.c_p ./ 2.0 .* a[xx] .^ 2.0 ./ b[xx] .+ (1.0 .- P.c_lambda) .* c[xx] .* b[xx]
    end
end

end

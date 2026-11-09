function fdepodrift!(x,a,b,c,d,P)
    N = length(x)
    y = 0.0
    z = 0.0
    for xx = 1:N
        y = a[xx]./b[xx] .-1.0
        if y > 0.0
            x[xx] = y .* c[xx]./P.c_p
            z = y .- P.c_p ./ 2.0 .* x[xx] ./ c[xx] .- (1.0 .- P.c_lambda) .* d[xx] .* c[xx] .* x[xx]
            if z < 0.0
                x[xx] = 0.0
            end
        else
            x[xx] = y .* c[xx]./P.c_n
            z = a[xx] .* x[xx] .- b[xx] .* (x[xx] .+ P.c_n ./ 2.0 .* x[xx] .^ 2.0 ./ c[xx] .+ (1.0 .- P.c_lambda) .* d[xx] .* c[xx] )
            if z < 0.0
                x[xx] = 0.0
            end
        end
    end
end


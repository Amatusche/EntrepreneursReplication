function depodrift!(x,a,b,c,χ0,χ1,χ2)

    for xx in eachindex(x)
        y = a[xx] ./ b[xx] .- 1.0
        if y > χ0
            x[xx] = ((y .- χ0) ./ (χ1 .* χ2)) .^ (1.0 ./ (χ2 .- 1.0)) .* c[xx]
        elseif y < - χ0
            x[xx] = -(-(y .+ χ0) ./ (χ1 .* χ2)) .^ (1.0 ./ (χ2 .- 1.0)) .* c[xx]
        else
            x[xx] = 0.0
        end
    end

end

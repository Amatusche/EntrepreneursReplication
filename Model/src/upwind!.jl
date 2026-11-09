function upwind!(d_B, dBF, dBB)

    for xx in eachindex(d_B)
        d_B[xx] = 0.0
        if dBF[xx] > 0.0
            d_B[xx] += dBF[xx]
        end
        if dBB[xx] < 0.0
            d_B[xx] += dBB[xx]
        end
    end
end

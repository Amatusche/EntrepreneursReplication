function sdebuild!(sde_B,d_B,a_grid_MAT_e_adj,χ0,χ1,χ2)
    for xx in eachindex(sde_B)

        if d_B[xx] > 0.0
            sde_B[xx] = -d_B[xx] .- χ0 .* d_B[xx] .- χ1 .* (d_B[xx] ./ a_grid_MAT_e_adj[xx]) .^ χ2 .* a_grid_MAT_e_adj[xx]
        else
            sde_B[xx] = -d_B[xx] .+ χ0 .* d_B[xx] .- χ1 .* (-d_B[xx] ./a_grid_MAT_e_adj[xx]) .^ χ2 .* a_grid_MAT_e_adj[xx]
        end
    end
end

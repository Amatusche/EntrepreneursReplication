"`change_distr(ge,share,factor,dim,grid)`

This function creates a distribution adjusted_ge over the `grid`.

Starting from the distribution ge, wealth of the richest `share` of the distribution along dimension dim increases by the factor `factor`."
function change_distr(ge,share,factor,dim,grid)
    # This function creates a distribution adjusted_ge over the grid in grid.
    # Starting from the distribution ge, wealth of the richest "share" of the distribution along dimension dim increases by the factor "factor".
    ndimensions = ndims(ge)
    dimvec = collect(1:ndimensions)
    deleteat!(dimvec, dim)
    N = size(ge)[dim]
    ge_marg = Array{Float64}(undef,N)
    ge_marg = reshape(sum(ge,dims=dimvec),N)
    prctl_index = findfirst(cumsum(ge_marg,dims=1).> 1.0.-share)

    adjusted_ge = similar(ge)
    selectdim(adjusted_ge,dim,1) .= selectdim(ge,dim,1)
    new_grid_point = Int64
    right_weight = Float64
    left_weight = Float64
    for ee = N:-1:max(prctl_index,2)
        new_grid_point = findfirst(min(grid[ee]*factor,maximum(grid)) .<= grid)
        right_weight = (min(grid[ee]*factor,maximum(grid)).- grid[new_grid_point.-1]) ./ (grid[new_grid_point].- grid[new_grid_point-1])
        left_weight = 1.0.-right_weight

        adjusted_ge[:,:,new_grid_point,:] .= adjusted_ge[:,:,new_grid_point,:] + right_weight*ge[:,:,ee,:];
        adjusted_ge[:,:,new_grid_point-1,:] .= adjusted_ge[:,:,new_grid_point-1,:] + left_weight*ge[:,:,ee,:];
    end
    for ee = prctl_index-1:-1:2

        adjusted_ge[:,:,ee,:] .= adjusted_ge[:,:,ee,:] .+ ge[:,:,ee,:]
    end
    adjusted_ge .= adjusted_ge./sum(adjusted_ge[:]);
end

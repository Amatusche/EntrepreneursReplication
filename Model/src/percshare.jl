"
`percshare(grid,pdf,p)`

Computes the share owned by the top `p`*100% for a discrete distribution
characterized by a `pmf` over an asset `grid`."
function percshare(grid, pmf, p)
    # Sort grid
    sortindex = sortperm(grid)
    grid .= grid[sortindex]
    pmf = pmf[sortindex]

    cmf = cumsum(pmf)
    index = findfirst(cmf .> (1.0 .- p))
    total = dot(grid, pmf)
    topp = dot(grid[index:end], pmf[index:end])
    return perc_share = topp ./ total

end

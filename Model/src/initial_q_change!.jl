function initial_q_change!(g, rel_change, grid, N, dim)

    a_grid_new = min.(rel_change .* grid, maximum(grid))
    g_new = zeros(size(g))

    intf = interpolate((grid,), 1:N, Gridded(Linear()))

    for xx = 1:N
        new_spot = intf(a_grid_new[xx])
        left_gridpoint = floor(Int64, new_spot)
        right_gridpoint = min(left_gridpoint+1, N)
        weight_on_leftpoint = 1.0 .- mod(new_spot, 1)
        weight_on_rightpoint = 1.0 .- weight_on_leftpoint
        if dim == 2
            g_new[:, left_gridpoint, :, :] = g_new[:, left_gridpoint, :, :] .+ weight_on_leftpoint*g[:, xx, :, :]
            g_new[:, right_gridpoint, :, :] = g_new[:, right_gridpoint, :, :] .+ weight_on_rightpoint*g[:, xx, :, :]
        elseif dim == 3
            g_new[:, :, left_gridpoint, :] = g_new[:, :, left_gridpoint, :] .+ weight_on_leftpoint*g[:, :, xx, :]
            g_new[:, :, right_gridpoint, :] = g_new[:, :, right_gridpoint, :] .+ weight_on_rightpoint*g[:, :, xx, :]
        end
    end

    g_new .= g_new ./ sum(g_new)*sum(g)

    # If rel_change is positive, some mass of the distribution will be pushed beyond the maximum gridpoint. Adjust for that
    old_sum = sum(grid[ii] * sum(selectdim(g, dim, ii)) for ii in 1:N)
    new_sum = sum(grid[ii] * sum(selectdim(g_new, dim, ii)) for ii in 1:N)
    gap = rel_change * old_sum - new_sum
    tol = 1e-12

    if gap > tol
        donor = selectdim(g_new, dim, N - 1)
        recipient = selectdim(g_new, dim, N)

        mass_to_move = gap / (grid[end] - grid[end-1])
        println("Adjusting for mass pushed beyond the maximum gridpoint in initial_q_change! Moving "* string(mass_to_move) * " units of mass from penultimate to last asset bin.")
        donor_mass = sum(donor)

        mass_to_move <= donor_mass || error("Insufficient mass in penultimate asset bin")

        fraction_to_move = mass_to_move / donor_mass

        recipient .+= fraction_to_move .* donor
        donor .*= 1.0 - fraction_to_move
    elseif gap < -tol
        error("Asset remapping unexpectedly increased the first moment")
    end

    g .= g_new

end

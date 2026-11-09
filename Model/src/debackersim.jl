function _fill_debacker_model!(model, NetProfits1, NetProfits2, Profits1, Profits2)
    NumberPercGroups = size(model, 1)
    N = length(NetProfits1)
    cutoffs1 = Vector{Float64}(undef, NumberPercGroups)
    cutoffs2 = Vector{Float64}(undef, NumberPercGroups)

    @inbounds for nn = 1:NumberPercGroups
        cutoff_idx = max(1, floor(Int, nn * N / NumberPercGroups))
        cutoffs1[nn] = Profits1[cutoff_idx]
        cutoffs2[nn] = Profits2[cutoff_idx]
    end

    fill!(model, 0.0)
    @inbounds for nn in eachindex(NetProfits1, NetProfits2)
        group1 = searchsortedfirst(cutoffs1, NetProfits1[nn])
        group2 = searchsortedfirst(cutoffs2, NetProfits2[nn])
        model[min(group1, NumberPercGroups), min(group2, NumberPercGroups)] += 1.0
    end

    @inbounds for nn = 1:NumberPercGroups
        rowsum = 0.0
        for mm = 1:NumberPercGroups
            rowsum += model[nn, mm]
        end
        if rowsum > 0.0
            for mm = 1:NumberPercGroups
                model[nn, mm] /= rowsum
            end
        else
            for mm = 1:NumberPercGroups
                model[nn, mm] = NaN
            end
        end
    end

    return model
end

" `debacker_sim(ge,Aentre,Ne,N1,N2,dt,TaxableBusinessIncomeE)`

Simulates a number N of households drawn from the distribution g for two years N2 times
and computes the matrix of transition probabilities between deciles of the
distribution of yearly taxable business income. `TaxableBusinessIncomeE` is a
quarterly flow; investment and changes in business equity are deliberately excluded.
"
function debackersim(ge, Aentre, Ne, N1, N2, dt, TaxableBusinessIncomeE)

    quarters_per_year = 4
    years_to_simulate = 2
    N = N1 * N2
    steps = Int(div(quarters_per_year * years_to_simulate, dt))
    steps_per_year = div(steps, years_to_simulate)

    initial_draws = Vector{Float64}(undef, N)
    @inbounds for nn = 1:N1
        initial_draws[nn] = rand()
    end
    @inbounds for rep = 2:N2
        copyto!(initial_draws, (rep - 1) * N1 + 1, initial_draws, 1, N1)
    end
    sort!(initial_draws)

    gcum = cumsum(vec(ge))
    states = Vector{Int32}(undef, N)
    next_states = similar(states)
    states_year1 = Vector{Int}(undef, N)
    states_year2 = Vector{Int}(undef, N)
    @inbounds for nn = 1:N
        states[nn] = Int32(searchsortedfirst(gcum, initial_draws[nn]))
    end

    println("   Simulating DeBacker et al. Transition Matrix: $N households for $steps time steps.")

    DeBackerModelFull = Array{Float64}(undef, 10, 10)
    DeBackerModelQuint = Array{Float64}(undef, 5, 5)

    AT_eplusI = I + dt .* Aentre'
    pstay = diag(AT_eplusI)
    colptr = AT_eplusI.colptr
    rowval = AT_eplusI.rowval
    nzval = AT_eplusI.nzval

    NetProfits1 = zeros(Float64, N)
    NetProfits2 = zeros(Float64, N)

    net_profit_base = TaxableBusinessIncomeE[:] .* dt

    for ss = 1:steps

        profits = ss <= steps_per_year ? NetProfits1 : NetProfits2
        @inbounds for nn = 1:N
            i = states[nn]
            first_idx = colptr[i]
            last_idx = colptr[i+1] - 1
            next_state = rowval[last_idx]
            draw = rand()
            if draw <= pstay[i] # stay in same state
                next_state = i
            else
                cumprob = pstay[i]
                for ptr in first_idx:last_idx
                    row = rowval[ptr]
                    if row == i # skip diagonal element (staying), because we already know that we leave the state
                        continue
                    end
                    cumprob += nzval[ptr]
                    if draw <= cumprob
                        next_state = rowval[ptr]
                        break
                    end
                end
            end
            next_states[nn] = Int32(next_state)
            profits[nn] += net_profit_base[i]
        end
        states, next_states = next_states, states
        if ss == steps_per_year
            copyto!(states_year1, states)
        elseif ss == steps
            copyto!(states_year2, states)
        end
    end

    # Distribution of simulated yearly profits (both years separately)
    copyto!(initial_draws, NetProfits1)
    sort!(initial_draws)
    Profits1 = initial_draws
    Profits2 = copy(NetProfits2)
    sort!(Profits2)

    _fill_debacker_model!(DeBackerModelQuint, NetProfits1, NetProfits2, Profits1, Profits2)
    _fill_debacker_model!(DeBackerModelFull, NetProfits1, NetProfits2, Profits1, Profits2)

    NetProfits = [NetProfits1; NetProfits2]
    states_current = [states_year1; states_year2]

    return DeBackerModelQuint, DeBackerModelFull, NetProfits, states_current

end

"""
Simulate two years of taxable pass-through business income in continuous time.

`TaxableBusinessIncomeE` is a quarterly flow. The simulation integrates this flow
along each state path without subtracting investment or adding capital gains.
"""
function debackersim2(ge, Aentre, N, TaxableBusinessIncomeE; rng=Random.default_rng())

    DeBackerModelFull = Matrix{Float64}(undef, 10, 10)
    DeBackerModelQuint = Matrix{Float64}(undef, 5, 5)
    NetProfits1 = zeros(Float64, N)
    NetProfits2 = zeros(Float64, N)
    states_year1 = Vector{Int}(undef, N)
    states_year2 = Vector{Int}(undef, N)

    profit_rate = vec(TaxableBusinessIncomeE)
    year_length = 4.0
    T = 8.0

    gcum = cumsum(vec(ge))
    Q = sparse(Aentre')
    jump_colptr, jump_rows, jump_cumrates, exit_rates = build_jump_cache(Q)

    @inbounds for nn = 1:N
        t = 0.0
        i = searchsortedfirst(gcum, rand(rng))
        stored_year1 = false

        while t < T
            exit_rate = exit_rates[i]

            if exit_rate <= 0.0
                # no more jumps; accrue profits until T
                accrue_profit!(NetProfits1, NetProfits2, nn, profit_rate[i], t, T)
                if !stored_year1 && t <= year_length
                    states_year1[nn] = i
                    stored_year1 = true
                end
                break
            end

            jump_time = t - log1p(-rand(rng)) / exit_rate
            t_next = min(jump_time, T)

            accrue_profit!(NetProfits1, NetProfits2, nn, profit_rate[i], t, t_next)
            if !stored_year1 && t <= year_length && year_length < t_next
                states_year1[nn] = i
                stored_year1 = true
            end

            if jump_time > T
                break
            end

            # choose destination conditional on a jump
            draw = rand(rng) * exit_rate
            next_i = jump_rows[jump_colptr[i+1]-1]

            @inbounds for ptr = jump_colptr[i]:(jump_colptr[i+1]-1)
                if draw <= jump_cumrates[ptr]
                    next_i = jump_rows[ptr]
                    break
                end
            end

            if !stored_year1 && jump_time == year_length
                states_year1[nn] = next_i
                stored_year1 = true
            end
            i = next_i
            t = jump_time
        end
        if !stored_year1
            states_year1[nn] = i
        end
        states_year2[nn] = i
    end

    # Distribution of simulated yearly profits (both years separately)
    Profits1 = copy(NetProfits1)
    Profits2 = copy(NetProfits2)
    sort!(Profits1)
    sort!(Profits2)

    _fill_debacker_model!(DeBackerModelQuint, NetProfits1, NetProfits2, Profits1, Profits2)
    _fill_debacker_model!(DeBackerModelFull, NetProfits1, NetProfits2, Profits1, Profits2)

    NetProfits = [NetProfits1; NetProfits2]
    states_current = [states_year1; states_year2]

    return DeBackerModelQuint, DeBackerModelFull, NetProfits, states_current

end

@inline function accrue_profit!(NetProfits1, NetProfits2, nn, rate, t0, t1)
    year_length = 4.0

    if t1 <= year_length
        NetProfits1[nn] += rate * (t1 - t0)
    elseif t0 >= year_length
        NetProfits2[nn] += rate * (t1 - t0)
    else
        NetProfits1[nn] += rate * (year_length - t0)
        NetProfits2[nn] += rate * (t1 - year_length)
    end
end

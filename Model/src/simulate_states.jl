import Random

"""
    simulate_states(initial_states, A, horizon; rng=Random.default_rng())

Simulate the continuous-time Markov chain with generator `A` for `horizon`
model-time units and return each household's final state. `A[i, j]` is the
transition rate from state `i` to state `j`, matching the convention used by
the HJB transition matrices in this model.

Passing the initial states explicitly makes it possible to start baseline and
counterfactual simulations from the same households. To reproduce random
draws, pass a seeded random-number generator through `rng`.
"""
function simulate_states(initial_states::AbstractVector{<:Integer}, A::AbstractMatrix,
    horizon::Real; rng::Random.AbstractRNG=Random.default_rng())

    states = Int.(initial_states)

    return _simulate_states!(states, A, horizon, rng)
end

"""
    simulate_states(initial_states, A_path, durations; rng=Random.default_rng())

Simulate a transition with a piecewise-constant generator. `A_path[n]` applies
for `durations[n]` model-time units. For an MIT transition, pass matching
slices of `A_t` (or `Aentre_t`) and `G.dt`.
"""
function simulate_states(initial_states::AbstractVector{<:Integer},
    A_path::AbstractVector{<:AbstractMatrix}, durations::AbstractVector{<:Real};
    rng::Random.AbstractRNG=Random.default_rng())

    length(A_path) == length(durations) ||
        throw(ArgumentError("A_path and durations must have the same length"))

    states = Int.(initial_states)
    for nn in eachindex(A_path, durations)
        _simulate_states!(states, A_path[nn], durations[nn], rng)
    end

    return states
end

function _simulate_states!(states, A, horizon, rng)
    # A is a backward generator (rows are origin states). Transposing it puts
    # the outgoing rates for one state in a sparse CSC column.
    Q = sparse(transpose(A))
    jump_colptr, jump_rows, jump_cumrates, exit_rates = build_jump_cache(Q)

    @inbounds for nn in eachindex(states)
        t = 0.0
        i = states[nn]

        while t < horizon
            exit_rate = exit_rates[i]
            exit_rate <= 0.0 && break

            jump_time = t + Random.randexp(rng) / exit_rate
            jump_time > horizon && break

            draw = rand(rng) * exit_rate
            first_ptr = jump_colptr[i]
            last_ptr = jump_colptr[i + 1] - 1
            next_i = jump_rows[last_ptr]
            for ptr = first_ptr:last_ptr
                if draw <= jump_cumrates[ptr]
                    next_i = jump_rows[ptr]
                    break
                end
            end

            i = next_i
            t = jump_time
        end

        states[nn] = i
    end

    return states
end

function build_jump_cache(Q)
    Ne = size(Q, 2)
    colptr = Q.colptr
    rowval = Q.rowval
    nzval = Q.nzval

    jump_colptr = Vector{Int}(undef, Ne + 1)
    jump_rows = Int[]
    jump_cumrates = Float64[]
    sizehint!(jump_rows, nnz(Q))
    sizehint!(jump_cumrates, nnz(Q))
    exit_rates = zeros(Float64, Ne)

    jump_colptr[1] = 1

    @inbounds for i = 1:Ne
        cumrate = 0.0
        for ptr = colptr[i]:(colptr[i+1]-1)
            j = rowval[ptr]
            rate = nzval[ptr]

            if j != i && rate > 0.0
                cumrate += rate
                push!(jump_rows, j)
                push!(jump_cumrates, cumrate)
            end
        end

        exit_rates[i] = cumrate
        jump_colptr[i+1] = length(jump_rows) + 1
    end

    return jump_colptr, jump_rows, jump_cumrates, exit_rates
end

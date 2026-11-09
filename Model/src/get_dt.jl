function get_dt(N_t,deltatransparam,deltatransmin,deltatransmax)

    x = Vector{Float64}(undef,N_t)
    x[1] = 0.0
    x[N_t] = 1.0
    for i = 2:N_t-1
        x[i] = (i-1)/(N_t-1)
    end

    z = x.^(1/deltatransparam)
    y = deltatransmin .+ (deltatransmax-deltatransmin).*z

end

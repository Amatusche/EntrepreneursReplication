function shock(shock1,shock_ac,N_t,dt)
    mp_shock_t = zeros(N_t);
    mp_shock_t[1] = shock1/4;
    for tt = 2:N_t
        mp_shock_t[tt] = mp_shock_t[tt-1].*shock_ac.^dt[tt-1];
    end
return mp_shock_t
end

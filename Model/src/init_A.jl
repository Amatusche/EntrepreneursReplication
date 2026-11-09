function init_A(N_b, N_a, N, Aexo)

    ## Initialize Aentre matrix
    Aendo = spdiagm(N,N,-N_b * N_a => zeros(N_b*N_a+1:N), -N_b => zeros(N_b+1:N), -1 => zeros(2:N),
        0 => zeros(N), 1 => zeros(1:N-1), N_b => zeros(1:N-N_b), N_b * N_a => zeros(1:N-N_b*N_a))
    i1, j1, vaexo1 = findnz(Aendo)
    i2, j2, vaexo2 = findnz(Aexo)
    iaexo = [i1; i2]
    jaexo = [j1; j2]
    vaexo = [vaexo1; vaexo2]
    A = sparse(iaexo, jaexo, vaexo)
    return A
end

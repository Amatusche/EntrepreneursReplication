function find_z_stat(N_j,beta_1,z_1,phi_1,lambda_1,delta_z1_B,delta_z1_F)

N = N_j;
N_k = 1;
N_ab = 1;
# This is the drift
Az1 = spzeros(N,N);
for k = 1:N_k
    for j = 1:N_j
        Az1[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1).*N_ab.*N_j.+j.*N_ab),((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab)] = sparse(I,N_ab,N_ab).*(- max(-beta_1*z_1[j],0)/delta_z1_F[j] .+ min(-beta_1*z_1[j],0)/delta_z1_B[j]);
        if j < N_j
            Az1[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+j*N_ab+1):((k-1)*N_ab*N_j+(j+1)*N_ab)] = sparse(I,N_ab,N_ab).*(max(-beta_1*z_1[j],0)/delta_z1_F[j]);
        end
        if j > 1
            Az1[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+(j-2)*N_ab+1):((k-1)*N_ab*N_j+(j-1)*N_ab)] = sparse(I,N_ab,N_ab).*(-min(-beta_1*z_1[j],0)/delta_z1_B[j]);
        end
    end
end

# This is the jump
Az3 = spzeros(N,N);
for k = 1:N_k
    for j = 1:N_j
        for jj = 1:N_j
            if j == jj
                Az3[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+(jj-1)*N_ab+1):((k-1)*N_ab*N_j+jj*N_ab)] = sparse(I,N_ab,N_ab).*(-lambda_1+lambda_1*phi_1[jj]);
            else
                Az3[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+(jj-1)*N_ab+1):((k-1)*N_ab*N_j+jj*N_ab)] = sparse(I,N_ab,N_ab).*(lambda_1*phi_1[jj]);
            end
        end
    end
end

Z = Az1+Az3;

# compute KFE
AT = Z';
b = zeros(N);
# need to fix one value, otherwise matrix is singular
iFix = 2;
b[iFix]=.1;
row = [zeros(iFix-1); 1; zeros(N-iFix)];
AT[iFix,:] = row;

 # Solve linear system
g_tilde = AT\b;
#rescale \tilde{g} so that it sums to 1
g_sum = sum(g_tilde);
g_tilde = g_tilde./g_sum;

g = reshape(g_tilde,N_j,N_k);

end

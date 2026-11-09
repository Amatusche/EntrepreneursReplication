function exostate_transit_mat(N_ab,N_j,N_k,beta_1,beta_2,lambda_1,lambda_2,phi_1,phi_2,z_1,z_2,delta_z1_F,delta_z1_B,delta_z2_F,delta_z2_B)

N = N_ab*N_j*N_k;
# This is the drift of z1
Az1 = spzeros(N,N);
for k = 1:N_k
    for j = 1:N_j
        Az1[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab)] = sparse(I,N_ab,N_ab).*(- max(-beta_1*z_1[j],0)/delta_z1_F[j] + min(-beta_1*z_1[j],0)/delta_z1_B[j]);
        if j < N_j
            Az1[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+j*N_ab+1):((k-1)*N_ab*N_j+(j+1)*N_ab)] = sparse(I,N_ab,N_ab).*(max(-beta_1*z_1[j],0)/delta_z1_F[j]);
        end
        if j > 1
            Az1[((k-1)*N_ab*N_j+(j-1)*N_ab+1):((k-1)*N_ab*N_j+j*N_ab),((k-1)*N_ab*N_j+(j-2)*N_ab+1):((k-1)*N_ab*N_j+(j-1)*N_ab)] = sparse(I,N_ab,N_ab).*(-min(-beta_1*z_1[j],0)/delta_z1_B[j]);
        end
    end
end

# This is the drift of z2
Az2 = spzeros(N,N);
for k = 1:N_k
    Az2[((k-1)*N_ab*N_j+1):(k*N_ab*N_j),((k-1)*N_ab*N_j+1):(k*N_ab*N_j)] = sparse(I,N_ab*N_j,N_ab*N_j).*(-max(-beta_2*z_2[k],0)/delta_z2_F[k] + min(-beta_2*z_2[k],0)/delta_z2_B[k]);
    if k < N_k
        Az2[((k-1)*N_ab*N_j+1):(k*N_ab*N_j),(k*N_ab*N_j+1):((k+1)*N_ab*N_j)] = sparse(I,N_ab*N_j,N_ab*N_j).*(max(-beta_2*z_2[k],0)/delta_z2_F[k]);
    end
    if k > 1
        Az2[((k-1)*N_ab*N_j+1):(k*N_ab*N_j),((k-2)*N_ab*N_j+1):((k-1)*N_ab*N_j)] = sparse(I,N_ab*N_j,N_ab*N_j).*(-min(-beta_2*z_2[k],0)/delta_z2_B[k]);
    end
end

# This is the jump of z1
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

# This is the jump of z2
Az4 = spzeros(N,N);
for k = 1:N_k
    for kk = 1:N_k
        if k == kk
            Az4[((k-1)*N_ab*N_j+1):(k*N_ab*N_j),((kk-1)*N_ab*N_j+1):(kk*N_ab*N_j) ]= sparse(I,N_ab*N_j,N_ab*N_j).*(-lambda_2+lambda_2*phi_2[kk]);
        else
            Az4[((k-1)*N_ab*N_j+1):(k*N_ab*N_j),((kk-1)*N_ab*N_j+1):(kk*N_ab*N_j)] = sparse(I,N_ab*N_j,N_ab*N_j).*(lambda_2*phi_2[kk]);
        end
    end
end

Az = Az1 + Az2 + Az3 + Az4;
end

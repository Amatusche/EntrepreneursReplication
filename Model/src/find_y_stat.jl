function find_y_stat(lambda_y_lh,lambda_y_hl)

N_y = 2;

# This is the jump
Ay = spzeros(N_y,N_y);

Ay[1,1] = -lambda_y_lh;
Ay[1,2] = lambda_y_lh;
Ay[2,2] = -lambda_y_hl;
Ay[2,1] = lambda_y_hl;

# compute KFE
AT = Ay';
b = zeros(N_y,1);
# need to fix one value, otherwise matrix is singular
iFix = 2;
b[iFix]=.001;
row = [zeros(iFix-1);1;zeros(N_y-iFix)];
AT[iFix,:] = row;

 # Solve linear system
g_tilde = AT\b;
#rescale \tilde{g} so that it sums to 1
g_sum = sum(g_tilde);
g_tilde = g_tilde./g_sum;

g = reshape(g_tilde,1,N_y);

end

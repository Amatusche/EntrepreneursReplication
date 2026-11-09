function exostate_transit_mat_entre(N_bae,N_y,lambda_y_lh,lambda_y_hl)

Ne = N_bae*N_y;

Ayll = spzeros(Ne,Ne);
Aylh = spzeros(Ne,Ne);
Ayhh = spzeros(Ne,Ne);
Ayhl = spzeros(Ne,Ne);

if N_y > 1
    Ayll[1:N_bae,1:N_bae] = sparse(I,N_bae,N_bae).*(-lambda_y_lh);
    Aylh[1:N_bae,(N_bae+1):(2*N_bae)] = sparse(I,N_bae,N_bae).*lambda_y_lh;
    Ayhh[(N_bae+1):(2*N_bae),(N_bae+1):(2*N_bae)] = sparse(I,N_bae,N_bae).*(-lambda_y_hl);
    Ayhl[(N_bae+1):(2*N_bae),1:N_bae] = sparse(I,N_bae,N_bae).*lambda_y_hl;
end

Ay = Ayll + Aylh + Ayhh + Ayhl;

end

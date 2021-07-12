


plot(repmat(x'*Diagrams,[N_sample 1])' + ((Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams)')

sigma_i = 0.5*sqrt(rho/N_sample*p(1));



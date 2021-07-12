
EVaR_range = [0:0.01:0.25];
%Optimal_mu_for_EVaR = [];
%Portfolios_nominal = [];
%Portfolios_robust  = [];

for i=17:length(EVaR_range)
    
    EVaR=EVaR_range(i);
    solver_recoded;
    %solver_file_nonrobust_otherway;
    
    Optimal_mu_for_EVaR(i,1) = value_robust;
    %Optimal_mu_for_EVaR(i,2) = value_nonrobust;
    
    %Portfolios_nominal = [Portfolios_nominal w_nonrobust];
    Portfolios_robust(:,i) = w_robust;
    
end
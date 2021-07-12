%%% EVaR's of robust and nonrobust portfolio for sampled probability
%%% distributions
EVaR_data=[];
%%% Mean returns of robust and nonrobust portfolio for sampled probability
%%% distributions
return_data=[];
%%% Phi divergence of the sampled probability distributions
dispersion_data=[];

% Main loop for the portfolio bootstrap test
for iterate_sample=1:1000
    portfolio_bootstrap_single_case;
    iterate_sample
end
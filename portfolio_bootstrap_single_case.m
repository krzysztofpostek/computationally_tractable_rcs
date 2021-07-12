%%%% This is the file samples a probability distribution and computes the
%%%% mean return and EVaR of the robust and the nonrobust portfolio

%%%% First we sample a probability distribution %%%%%%%%%%%%%
found_probability=0;
std_multiplier=1;

while(found_probability ==0)
    probability_trial=1/N+randn(N-1,1)*std_multiplier*sqrt(rho_Q/N^2);
    probability_trial=[probability_trial; 1-sum(probability_trial)];
    if(prod(probability_trial>=0))
        found_probability=1;
        probabilities=probability_trial;
    end
end

%%% We compute the value of the divergence of the sampled probability
%%% measure relative to the uniform one
phi_divergence_value=sum(((probabilities-ones(N,1)/N).^2)./(ones(N,1)/N));

%%% Now we find the EVaR value of the robust portfolio
cvx_begin quiet
    variable lambda nonnegative;
    variables eta v(N) z;
    minimize(z)
    subject to
        % Main constraint
        eta+rho*lambda - lambda + sum(probabilities.*v) <= z;
        
        for n=1:N
            {-Y(:,n)'*w_robust - eta,lambda,v(n) } == exponential;
        end    
cvx_end

EVaR_robust=z;
%%% Computing the mean return of the robust portfolio
return_robust=w_robust'*Y*probabilities;

%%% Now we find the EVaR value of the nonrobust portfolio

cvx_begin quiet
    variable lambda nonnegative;
    variables eta v(N) z;
    minimize(z)
    subject to
        % Main constraint
        eta+rho*lambda - lambda + sum(probabilities.*v) <= z;
        
        for n=1:N
            {-Y(:,n)'*w_nonrobust - eta,lambda,v(n) } == exponential;
        end
        
cvx_end

EVaR_nonrobust=z;
%%% Computing the mean return of the nonrobust portfolio
return_nonrobust=w_nonrobust'*Y*probabilities;

%%% Adding the data
EVaR_data=[EVaR_data; EVaR_robust EVaR_nonrobust];
return_data=[return_data; return_robust return_nonrobust];
dispersion_data=[dispersion_data; phi_divergence_value];

save('Bootstrap_results_final_version.mat');
% Bootstrappin data generation

% Bootstrap test

N_sample_bootstrap_outer = 100;
N_sample_bootstrap_inner = 1000;

Sample_indices = zeros(N_sample_bootstrap_outer*N_sample_bootstrap_inner,1);

counter = 0;

rho_for_bootstrap = 0.1;

sigma_i = 0.5*rho_for_bootstrap/N_sample;
p_i = p(1);

probability_thing = 0;

while(counter < N_sample_bootstrap_outer)
    
    probability_sampled = p(1) + randn(N_sample-1,1)*sigma_i;
    probability_sampled = [probability_sampled; 1-sum(probability_sampled)];
    
    if(prod(probability_sampled >= 0) == 1)
        
        probability_thing = probability_thing + (sum(abs(probability_sampled - p)) <= rho_for_bootstrap);
    
        for iterate_sample = 1:N_sample_bootstrap_inner
            
            Sample_indices(counter*N_sample_bootstrap_inner +iterate_sample) = randsample([1:N_sample]',1,true,probability_sampled);
            
        end

        counter = counter + 1;
    
    end
end
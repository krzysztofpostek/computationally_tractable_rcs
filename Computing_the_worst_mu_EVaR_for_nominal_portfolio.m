%%%%%%%%%%%%%%%%%%% Determining the worst case EVaR %%%%%%%%%%%%

Nominal_portfolio_wc_EVaR_corresponding_return = [];
Nominal_portfolio_wc_return_corresponding_EVaR = [];

for i = 1:length(EVaR_range)
    
    w_nonrobust = Portfolios_nominal(:,i);
    
    cvx_begin
    
        variable q(N) nonnegative
        variable probsy(N) nonnegative
        variable z;
        maximize(z)
        subject to
        
            % Main constraint
            z <= sum(probsy.*(-(Y')*w_nonrobust));

            sum(rel_entr(probsy,q)) <= rho;
            sum(((q-r').^2)./(r')) <= rho_Q;

            q >= 0;
            probsy >= 0;

            sum(q) == 1;
            sum(probsy) == 1;

    cvx_end
    
    Nominal_portfolio_wc_EVaR_corresponding_return = [Nominal_portfolio_wc_EVaR_corresponding_return; probsy'*(-(Y')*w_nonrobust) q'*((Y')*w_nonrobust)];
   
    
end

for i=1:length(EVaR_range)

    w_nonrobust= Portfolios_nominal(:,i);

    cvx_begin
        variables probsy(N) z;
        minimize(z)
        subject to
            % Main constraint
            sum(probsy.*(Y'*w_nonrobust))<=z;

            sum(((probsy-r').^2)./(r')) <= rho_Q;

            probsy >= 0;

            ones(1,N)*probsy == 1;

    cvx_end
    
    cvx_begin
        variables q(N) z;
        maximize(z)
        subject to
            % Main constraint
            z <= sum(q.*(-(Y')*w_nonrobust));

            sum(rel_entr(q,probsy)) <= rho;

            q>=0;

            ones(1,N)*q==1;

    cvx_end

    Nominal_portfolio_wc_return_corresponding_EVaR = [Nominal_portfolio_wc_return_corresponding_EVaR; cvx_optval probsy'*(Y'*w_nonrobust)];
    
end

Robust_portfolio_wc_EVaR_corresponding_return = [];
Robust_portfolio_wc_return_corresponding_EVaR = [];

for i=1:length(EVaR_range)
    
    w_robust = Portfolios_robust(:,i);
    
    cvx_begin
        variable q(N) nonnegative 
        variable probsy(N) nonnegative
        variable z;
        maximize(z)
        subject to
        
            % Main constraint
            z <= sum(probsy.*(-(Y')*w_robust));

            sum(rel_entr(probsy,q)) <= rho;
            sum(((q-r').^2)./(r')) <= rho_Q;

            sum(q) ==1;
            sum(probsy) ==1;

    cvx_end
    
    Robust_portfolio_wc_EVaR_corresponding_return = [Robust_portfolio_wc_EVaR_corresponding_return; probsy'*(-(Y')*w_robust) q'*((Y')*w_robust)];
   
    
end

for i=1:length(EVaR_range)

    w_robust= Portfolios_robust(:,i);

    cvx_begin
        variables probsy(N) z;
        minimize(z)
        subject to
            % Main constraint
            sum(probsy.*(Y'*w_robust))<=z;

            sum(((probsy-r').^2)./(r')) <= rho_Q;

            probsy>=0;

            ones(1,N)*probsy==1;

    cvx_end
    
    cvx_begin
        variables q(N) z;
        maximize(z)
        subject to
            % Main constraint
            z <= sum(q.*(-(Y')*w_robust));

            sum(rel_entr(q,probsy)) <= rho;

            q>=0;

            ones(1,N)*q==1;

    cvx_end

    Robust_portfolio_wc_return_corresponding_EVaR = [Robust_portfolio_wc_return_corresponding_EVaR; cvx_optval probsy'*(Y'*w_robust) ];
    
end


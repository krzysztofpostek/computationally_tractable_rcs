%%% This file, for a given EVaR parameter value, computes the maximum mean
%%% return in the model with uncertainty in q

cvx_begin

    variable u(7) nonnegative;
    variable w(M) nonnegative;
    variables v(N) v1(N) v2(N) v3(2*N) v4(N) v5(N) v6(N) v7(N) z1(N) z2(N) z3(N) eta muvar;
    
    maximize(muvar)
    
    subject to
    
        u(1) - u(2) + u(3)*rho + u(4) - u(5) + u(6)*rho_Q + r*z1 <= EVaR;

        for n=1:N
            
            % Constraint on the EVaR
            max( -u(6) , v6(n) + 0.25 *quad_over_lin(v6(n),u(6)) ) <= z1(n);
            
            {v3(n),u(3),z3(n)} == exponential(1);
            
            % Constraint on the mean return
            max(-u(7) , v7(n)-eta + 0.25*quad_over_lin(v7(n)-eta,u(7)) ) <= z2(n);
        
        end

        v3(N+1:2*N) - u(3)*ones(N,1) + z3 <= 0;

        v1 <= u(1)*ones(N,1);
        v2 <= -u(2)*ones(N,1);

        v4 <= u(4)*ones(N,1);
        v5 <= -u(5)*ones(N,1);
        
        v3(1:N) + v4 + v5 + v6 == v;
        v1 + v2 + v3(N+1:2*N) == 0;

        -Y'*w <= v;
        -Y'*w <= v7;
        
        eta+u(7)*rho_Q +sum(r'.*z2) <= -muvar;

        sum(w) == 1;
    
cvx_end 

value_robust=cvx_optval;
w_robust=w;

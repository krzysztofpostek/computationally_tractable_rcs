cvx_begin

    variable u(7,1) nonnegative
    variable w(M) nonnegative
    variables z1(N) z2(N) z3(N)
    variables v(N) v1(N) v2(2*N) v3(2*N) v4(2*N) v5(2*N) v6(2*N) v7(2*N) eta mean_return

    maximize(mean_return)
    subject to 
    
        v2(1:N) <= u(2);
        v3(1:N) <= -u(3);
        
        v5(N+1:2*N) <= u(5);
        v6(N+1:2*N) <= -u(6);
        
        v2 + v3 + v4 + v5 + v6 + v7 == [v ; zeros(N,1)];
        
        v2(N+1:2*N) <= 0 ;
        v3(N+1:2*N) <= 0 ;
        
        v5(1:N) <= 0;
        v6(1:N) <= 0;
        v7(1:N) <= 0;
        
        -Y'*w <= v;
        -Y'*w <= v1;
        
        eta + u(1)*rho_Q + r*z1 <= -mean_return;
        
        u(2) - u(3) + u(4)*rho + u(5) - u(6) + u(7)*rho_Q + r*z2 <= EVaR;
        
        for n = 1:N
            
            max(-u(7), v7(N+n) + 0.25*quad_over_lin(v7(N+n),u(7)) ) <= z2(n);
            
            max(-u(1), v1(n) - eta + 0.25*quad_over_lin( v1(n) - eta , u(1)) ) <= z1(n);
            
            v4(N+n) + z3(n) - u(4) <= 0;
            {v4(n), u(4), z3(n)} == exponential;
            
        end
        
        sum(w) == 1;

cvx_end

value_robust = cvx_optval;
w_robust = w;
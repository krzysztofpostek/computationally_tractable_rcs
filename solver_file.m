% This file, for a given value of minimum mean return mu, computes the
% robust portfolio with maximum mean return

tic
cvx_begin
    variable u(7) nonnegative;
    variable w(M) nonnegative;
    variable q(N) nonnegative;
    variables v(N) v1(2*N) v2(2*N) v3(2*N) v4(2*N) v5(2*N) v6(2*N) v7(N);
    variables z1(N) z2(N) z eta;
    
    minimize(z)
    
    subject to
    u(1)-u(2)+u(3)*rho+u(4)-u(5)+u(6)*rho_Q+ sum(r'.*z1) <= z;
    
    for n=1:N
        max(-u(6),v6(N+n) + 0.25 *quad_over_lin(v6(N+n),u(6))) <=z1(n);
    end
    
    for n=1:N
        {v3(n),u(3),q(n)} == exponential;
    end
    
    for n=1:N
        v3(N+n) - u(3) + q(n) <=0;
    end
    
    v1(1:N) <= u(1)*ones(N,1);
    v2(1:N) <= -u(2)*ones(N,1);
    v4(N+1:2*N) <= u(4)*ones(N,1);
    v5(N+1:2*N) <= -u(5)*ones(N,1);
    
    v1(N+1:2*N)==zeros(N,1);
    v2(N+1:2*N)==zeros(N,1);
    v4(1:N)==zeros(N,1);
    v5(1:N)==zeros(N,1);
    v6(1:N)==zeros(N,1);
    
    v1+v2+v3+v4+v5+v6==AC'*v;
    
    for n=1:N
        -Y(:,n)'*w <= v(n);
        -Y(:,n)'*w == v7(n);
    end
    
    eta+u(7)*rho_Q+sum(r'.*z2) <= -mu;
    
    ones(1,M)*w==1;
    
    for n=1:N
        max(-u(7),v7(n)-eta + 0.25*quad_over_lin(v7(n)-eta,u(7)))<=z2(n);
    end
cvx_end

time_robust=toc;
value_robust=cvx_optval;
% This file, for a given value of minimum mean return mu, computes the
% nonrobust portfolio with maximum mean return
tic
cvx_begin
    variable u nonnegative; 
    variable w(M) nonnegative;
    variables z eta z1(N) z2(N) v1(N);
    minimize(z)
    subject to
    %Constraint 1
        eta + u*rho - u +sum(r'.*z1) <=z;
    %Constraint 2
        for n=1:N
             {v1(n)-eta,u,z1(n)} == exponential;
        end
    %Constraint 3
        for n=1:N
            -Y(:,n)'*w <= v1(n)
        end
    %Constraint 4
        sum(r'.*z2) <= -mu;
    %Constraint 5
        for n=1:N
            -Y(:,n)'*w <= z2(n)
        end
    %Constraint 6
        ones(M,1)'*w==1;
cvx_end

time_nonrobust=toc;
value_robust=cvx_optval;
w_robust=w;
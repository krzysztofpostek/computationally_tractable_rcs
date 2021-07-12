%%%%% This file for a given paramter EVaR computes the maximum mean return
%%%%% in a model without uncertainty
cvx_begin
    variable u nonnegative; 
    variable w(M) nonnegative;
    variables eta z1(N) v1(N) muvar;
    maximize(r*Y'*w)
    subject to
    %Constraint 1
        eta + u*rho - u +sum(r'.*z1) <=EVaR;
        
    %Constraint 2
        for n=1:N
             {v1(n)-eta,u,z1(n)} == exponential;
        end
        
    %Constraint 3
        -Y'*w <= v1;
            
    %Constraint 6
        sum(w) == 1;
        
cvx_end

value_nonrobust=cvx_optval;
w_nonrobust=w;
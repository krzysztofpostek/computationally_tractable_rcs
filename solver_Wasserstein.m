
D = zeros(4*N_sample + 3 , N_sample^2 + N_sample);
d = zeros(4*N_sample + 3,1 );

D(1,1:N_sample) = 1;
d(1) = 1;
D(2,1:N_sample) = -1;
d(2) = -1;

for i=1:N_sample
    for n = 1:N_sample
        D(3,N_sample*i + n) = norm(Error_matrix(:,i) - Error_matrix(:,i))^2;
    end
end

d(3) = rho;


for i=1:N_sample
    for n = 1:N_sample
        D(3+n,n) = -1;
        D(3+n,N_sample*n+i) = 1;
        
        D(3+N_sample+n,n) = 1;
        D(3+N_sample+n,N_sample*n+i) = -1;
        
        D(3+2*N_sample+n,N_sample*i+n) = 1;
        d(3+2*N_sample+n) = p(n);
        
        D(3+3*N_sample+n,N_sample*i+n) = -1;
        d(3+3*N_sample+n) = -p(n);
    end
end


cvx_begin
    variable tau
    variables kappa_left(N_grid_points,1) kappa_right(N_grid_points,1)
    variables v_left(N_sample,N_grid_points) v_right(N_sample,N_grid_points)
    variable u_left(4*N_sample+3,N_grid_points) nonnegative
    variable u_right(4*N_sample+3,N_grid_points) nonnegative
    variable x(N_antennas)
    minimize(tau)
    subject to
        
        for iterate_grid = min(Indices_class_1):max(Indices_class_1)
            
           -kappa_right(iterate_grid) + u_right(:,iterate_grid)'*d <= 0;
           
           [v_right(:,iterate_grid); zeros(N_sample^2,1)] <= D'*u_right(:,iterate_grid);
           
           v_right(:,iterate_grid) >=  1/alpha*max(0 , - tau + kappa_right(iterate_grid) + x'*Diagrams(:,iterate_grid) + (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );
           
           
           
           -kappa_left(iterate_grid) + u_left(:,iterate_grid)'*d <= 0;
           
           [v_left(:,iterate_grid); zeros(N_sample^2,1)] <= D'*u_left(:,iterate_grid);
           
           v_left(:,iterate_grid) >=  1/alpha*max(0 , - tau + kappa_left(iterate_grid) - x'*Diagrams(:,iterate_grid) - (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );

        end
        
        for iterate_grid = min(Indices_class_2):max(Indices_class_2)
            
           -kappa_right(iterate_grid) + u_right(:,iterate_grid)'*d <= 0;
           
           [v_right(:,iterate_grid); zeros(N_sample^2,1)] <= D'*u_right(:,iterate_grid);
           
           v_right(:,iterate_grid) >=  1/alpha*max(0 , - 1 + kappa_right(iterate_grid) + x'*Diagrams(:,iterate_grid) + (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );
           
           
           
           -kappa_left(iterate_grid) + u_left(:,iterate_grid)'*d <= 0;
           
           [v_left(:,iterate_grid); zeros(N_sample^2,1)] <= D'*u_left(:,iterate_grid);
           
           v_left(:,iterate_grid) >=  1/alpha*max(0 , - 1 + kappa_left(iterate_grid) - x'*Diagrams(:,iterate_grid) - (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );

        end
        
        for iterate_grid = min(Indices_class_3):max(Indices_class_3)
            
           -kappa_right(iterate_grid) + u_right(:,iterate_grid)'*d <= 0;
           
           [v_right(:,iterate_grid); zeros(N_sample^2,1)] <= D'*u_right(:,iterate_grid);
           
           v_right(:,iterate_grid) >=  1/alpha*max(0 , -1 + kappa_right(iterate_grid) + x'*Diagrams(:,iterate_grid) + (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );
           
           
           
           -kappa_left(iterate_grid) + u_left(:,iterate_grid)'*d <= 0;
           
           [v_left(:,iterate_grid); zeros(N_sample^2,1)] <= D'*u_left(:,iterate_grid);
           
           v_left(:,iterate_grid) >=  1/alpha*max(0 , 0.9 + kappa_left(iterate_grid) - x'*Diagrams(:,iterate_grid) - (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );

        end

cvx_end
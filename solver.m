% This file solves the antenna problem with correlated antenna errors

cvx_begin
    variable tau
    variables kappa_left(N_grid_points,1) kappa_right(N_grid_points,1)
    variables v_left(N_sample,N_grid_points) v_right(N_sample,N_grid_points)
    variable u_left(N_grid_points,1) nonnegative
    variable u_right(N_grid_points,1) nonnegative
    variable eta_left(N_grid_points,1)
    variable eta_right(N_grid_points,1)
    variable x(N_antennas)
    minimize(tau)
    subject to
        
        for iterate_grid = min(Indices_class_1):max(Indices_class_1)
            
           -kappa_right(iterate_grid) + eta_right(iterate_grid) + u_right(iterate_grid)*rho + p'*max(-u_right(iterate_grid) , v_right(:,iterate_grid) - eta_right(iterate_grid)  ) <= 0;
           
           v_right(:,iterate_grid) - eta_right(iterate_grid) <= u_right(iterate_grid);
           
           v_right(:,iterate_grid) >=  1/alpha*max(0 , - tau + kappa_right(iterate_grid) + x'*Diagrams(:,iterate_grid) + (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );
           
           
           
           -kappa_left(iterate_grid) + eta_left(iterate_grid) + u_left(iterate_grid)*rho + p'*max(-u_left(iterate_grid) , v_left(:,iterate_grid) - eta_left(iterate_grid)  ) <= 0;
           
           v_left(:,iterate_grid) - eta_left(iterate_grid) <= u_left(iterate_grid);
           
           v_left(:,iterate_grid) >=  1/alpha*max(0 , - tau + kappa_left(iterate_grid) - x'*Diagrams(:,iterate_grid) - (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );

        end
        
        for iterate_grid = min(Indices_class_2):max(Indices_class_2)
            
           -kappa_right(iterate_grid) + eta_right(iterate_grid) + u_right(iterate_grid)*rho + p'*max(-u_right(iterate_grid) , v_right(:,iterate_grid) - eta_right(iterate_grid)  ) <= 0;
           
           v_right(:,iterate_grid) - eta_right(iterate_grid) <= u_right(iterate_grid);
           
           v_right(:,iterate_grid) >=  1/alpha*max(0 , - 1 + kappa_right(iterate_grid) + x'*Diagrams(:,iterate_grid) + (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );
           
           
           
           -kappa_left(iterate_grid) + eta_left(iterate_grid) + u_left(iterate_grid)*rho + p'*max(-u_left(iterate_grid) , v_left(:,iterate_grid) - eta_left(iterate_grid)  ) <= 0;
           
           v_left(:,iterate_grid) - eta_left(iterate_grid) <= u_left(iterate_grid);
           
           v_left(:,iterate_grid) >=  1/alpha*max(0 , - 1 + kappa_left(iterate_grid) - x'*Diagrams(:,iterate_grid) - (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );

        end
        
        for iterate_grid = min(Indices_class_3):max(Indices_class_3)
            
           -kappa_right(iterate_grid) + eta_right(iterate_grid) + u_right(iterate_grid)*rho + p'*max(-u_right(iterate_grid) , v_right(:,iterate_grid) - eta_right(iterate_grid)  ) <= 0;
           
           v_right(:,iterate_grid) - eta_right(iterate_grid) <= u_right(iterate_grid);
           
           v_right(:,iterate_grid) >=  1/alpha*max(0 , -1 + kappa_right(iterate_grid) + x'*Diagrams(:,iterate_grid) + (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );
           
           
           
           -kappa_left(iterate_grid) + eta_left(iterate_grid) + u_left(iterate_grid)*rho + p'*max(-u_left(iterate_grid) , v_left(:,iterate_grid) - eta_left(iterate_grid)  ) <= 0;
           
           v_left(:,iterate_grid) - eta_left(iterate_grid) <= u_left(iterate_grid);
           
           v_left(:,iterate_grid) >=  1/alpha*max(0 , 0.9 + kappa_left(iterate_grid) - x'*Diagrams(:,iterate_grid) - (Error_matrix.*repmat(x,[1 N_sample]))'*Diagrams(:,iterate_grid) );

        end

cvx_end

% This file solves the antenna problem - nominal version

cvx_begin
    variable tau
    variables kappa_left(N_grid_points,1) kappa_right(N_grid_points,1)
    variable x(N_antennas)
    minimize(tau)
    subject to
        
        for iterate_grid = min(Indices_class_1):max(Indices_class_1)
            
           x'*Diagrams(:,iterate_grid)  <= tau;
           
           x'*Diagrams(:,iterate_grid)  >= -tau ;

        end
        
        for iterate_grid = min(Indices_class_2):max(Indices_class_2)
            
           x'*Diagrams(:,iterate_grid) <= 1;
           
           x'*Diagrams(:,iterate_grid) >= -1;

        end
        
        for iterate_grid = min(Indices_class_3):max(Indices_class_3)
            
           x'*Diagrams(:,iterate_grid)  <= 1;
           
           x'*Diagrams(:,iterate_grid)  >= 0.9;

        end

cvx_end

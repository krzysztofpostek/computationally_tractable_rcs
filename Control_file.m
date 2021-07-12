% Control file

Solutions = [];
Objective = [];
Probability_of_violation = [];
Probability_of_violation_normal = [];

% First we do the solution with no error

iterate_solution = 1;

Data_setting;
%solver_nominal;

%Solutions(:,iterate_solution) = x;
%Objective(iterate_solution) = cvx_optval;

% Then solutions with error

Bootstrap_normal_generation;
Bootstrap_generation;

for rho = 0:0.01:0.1
    
    iterate_solution = iterate_solution + 1;
     
    solver;
     
    Solutions(:,iterate_solution) = x ;
    Objective(iterate_solution) = cvx_optval;
    
    Bootstrappin;
    Probability_of_violation(iterate_solution) = violated;
     
    Bootstrappin_normal;
    Probability_of_violation_normal(iterate_solution) = violated;
    
end


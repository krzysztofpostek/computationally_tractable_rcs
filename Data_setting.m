
% Specify the density of grid on [0,90]
N_grid_points = 50;
Angle_grid = linspace(0,pi/2,N_grid_points);
Degrees_angle_grid = linspace(0,90,N_grid_points);
N_antennas = 40;
lambda = 0.25;

cvx_solver Gurobi

% Separation of the indices corresponding to different sectors of the
% [0,90] interval

Indices_class_1 = sum(Angle_grid < 70/180*pi);
Indices_class_1 = [1:Indices_class_1]';

Indices_class_2 = sum(Angle_grid < 77/180*pi);
Indices_class_2 = [length(Indices_class_1)+1:Indices_class_2]';

Indices_class_3 = N_grid_points;
Indices_class_3 = [length(Indices_class_1)+length(Indices_class_2)+1:N_grid_points]';

% Creating the vectors of diagrams of the k-th antenna for each of the
% point in the grid

Diagrams = zeros(N_grid_points,40);

for iterate_antenna = 1:40
    Energy = @(fi) 0.5*cos(2*pi*iterate_antenna/lambda/40*cos(Angle_grid)*cos(fi));
    
    Diagrams(:,iterate_antenna) = integral(Energy,0,2*pi,'ArrayValued',true)';
end

Diagrams = Diagrams';


N_sample = 100;
std_deviation_common = 0.005;
std_deviation_intrinsic = 0.005;

p = ones(N_sample,1)/N_sample;

Error_matrix = std_deviation_common*repmat(randn(1,N_sample),[N_antennas 1]) + std_deviation_intrinsic*randn(N_antennas,N_sample);

alpha = 0.1;

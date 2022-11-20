% Get the coordinates of the airfoil - airfoiltools.com - http://airfoiltools.com/airfoil/details?airfoil=goe802-il
airfoil_name = 'GOE802';
airfoil_name = 'NACA0009';
airfoil_coordinates_file = sprintf('%s_DiscretePoints.csv', airfoil_name);
GOE802_DiscretePoints = readtable(airfoil_coordinates_file, 'NumHeaderLines', 1);
GOE802_DiscretePoints = GOE802_DiscretePoints{:, :};

% Put coordinates in a single vector starting from trailing edge and ending there - going lower -> upper edge
x_all_original = GOE802_DiscretePoints(:, 1);
y_all_original = GOE802_DiscretePoints(:, 2);

% Increasing the fineness of the discretization
N_airfoil = 1000;

% Get the finer distribution and overwrite the original coordinates
% Using interparc function by John D Errico from File Exchange - https://in.mathworks.com/matlabcentral/fileexchange/34874-interparc
Airfoil_points_fine = interparc(N_airfoil, x_all_original, y_all_original, 'spline');
x_all_original = Airfoil_points_fine(:, 1);
y_all_original = Airfoil_points_fine(:, 2);

FineDiscretization_Coordinates = [x_all_original, y_all_original];

% Save airfoil points to file
filename = sprintf('%s_DiscretePoints_%d.mat', airfoil_name, N_airfoil);
save(filename, 'FineDiscretization_Coordinates');



















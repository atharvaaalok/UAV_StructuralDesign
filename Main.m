%% AS5100 Mini Project - Structural Design
% Authors: Team SARAS
% Last Modified: 28/ 08/ 2022


%% INTIAL SETTINGS
clc; clear; close all;
PS = PLOT_STANDARDS();


%% AIRCRAFT PARAMETERS, GENERAL CONSTANTS

global c_root c_tip AR b S lambda W g
c_root = 0.3731;                                % Root chord
c_tip = 0.1679;                                 % Tip chord
thickness_ratio = 9.8 / 100;                    % thickness ratio of airfoil GOE802
AR = 7.2;                                       % Aspect ratio of wing
b = 1.9475;                                     % Wing span
S = 2 * ((1/2) * (b/2) * (c_root + c_tip));     % Area of Wing
lambda = 0.45;                                  % Taper ratio of wing
W = 10.47;                                      % Total weight of UAV
g = 9.81;


%% AIRFOIL CURVE AND INCREASING DATA POINTS

load('GOE802_DiscretePoints_1000.mat');
% load('NACA0009_DiscretePoints_1000.mat');

% Get the fine discretization of the airfoil
x_all_original = FineDiscretization_Coordinates(:, 1);
y_all_original = FineDiscretization_Coordinates(:, 2);

% theta = [360; 330; 270; 225; 180; 120; 90; 45; 0];
% theta = theta * (pi / 180);
% x_all_original = cos(theta);
% y_all_original = sin(theta);

% Plot the airfoil
figure(1);
hold on; grid on; axis equal
plot(x_all_original, y_all_original);

% Increasing the fineness of the discretization
N_airfoil = length(x_all_original);

% interparc creates panels of equal length
airfoil_discretization_panel_length = sqrt( (x_all_original(2) - x_all_original(1))^2 + (y_all_original(2) - y_all_original(1))^2 );

% Put centroid at origin
centroid_x = sum(x_all_original) / length(x_all_original);
centroid_y = sum(y_all_original) / length(y_all_original);
x_all_original = x_all_original - centroid_x;
y_all_original = y_all_original - centroid_y;

% Find thickness of the airfoil and the location of max thickness. Spar will be placed here.
[x_min, x_min_Idx] = min(x_all_original);
x_lower = x_all_original(1: x_min_Idx);
y_lower = y_all_original(1: x_min_Idx);
x_upper = x_all_original(x_min_Idx + 1: end);
y_upper = y_all_original(x_min_Idx + 1: end);

N_compare = min([length(x_lower), length(x_upper)]);

airfoil_thickness = 0;
max_thickness_Idx = 1;
for i_compare = 1: N_compare
    x_lower_compare = x_lower(i_compare);
    [minVal, minVal_Idx] = min(abs(x_upper - x_lower_compare));
    x_upper_compare = x_upper(minVal_Idx);
    y_lower_compare = y_lower(i_compare);
    y_upper_compare = y_upper(minVal_Idx);

    temp = y_upper_compare - y_lower_compare;
    if temp > airfoil_thickness
        airfoil_thickness = temp;
        spar_Idx_1 = i_compare;
        spar_Idx_2 = length(x_lower) + minVal_Idx;
    end

end

x_max_thickness_lower = x_all_original(spar_Idx_1);
y_max_thickness_lower = y_all_original(spar_Idx_1);
x_max_thickness_upper = x_all_original(spar_Idx_2);
y_max_thickness_upper = y_all_original(spar_Idx_2);

x_spar_original = [x_max_thickness_lower; x_max_thickness_upper];
y_spar_original = [y_max_thickness_lower; y_max_thickness_upper];

% Break the spar into small panels similar to airfoil.
% Increasing the fineness of the discretization - keep panel length same
N_spar = round(sqrt( (x_max_thickness_upper - x_max_thickness_lower)^2 + (y_max_thickness_upper - y_max_thickness_lower)^2 ) / airfoil_discretization_panel_length);

% Get the finer distribution and overwrite the original coordinates
Spar_points_fine = interparc(N_spar, x_spar_original, y_spar_original, 'spline');
x_spar_original = Spar_points_fine(:, 1);
y_spar_original = Spar_points_fine(:, 2);

% % Rotate the airfoil points by an angle
% theta = pi/2;
% Rotation_Matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
% A = [x_all_original, y_all_original];
% B = [x_spar_original, y_spar_original];
% 
% A = A * Rotation_Matrix;
% B = B * Rotation_Matrix;
% 
% x_all_original = A(:, 1);
% y_all_original = A(:, 2);
% x_spar_original = B(:, 1);
% y_spar_original = B(:, 2);

% Getting coordinates at the root
chord_original = max(x_all_original) - min(x_all_original);
x_root = x_all_original * (c_root / chord_original);
y_root = y_all_original * (c_root / chord_original);

x_spar_root = x_spar_original * (c_root / chord_original);
y_spar_root = y_spar_original * (c_root / chord_original);

figure(2);
hold on; grid on; axis equal
plot(x_root, y_root);
plot(x_spar_root, y_spar_root);
plot(0, 0, 'og');


%% SCHRENK'S METHOD FOR LIFT DISTRIBUTION - SHEAR FORCE AND BENDING MOMENT DIAGRAMS

% Define number of points along wing span
N_z = 1000;
z = linspace(0, b/2, N_z);
dz = z(2) - z(1);

% Get chord distribution of wing
c_z = c_root * Length_ScalingFactor_at_z(z, b, lambda);

% Get chord distribution of an elliptic wing with same Area and Span
c_elliptic_z = ((4 * S) / (pi * b)) * sqrt(1-(2*z/b).^2);

% Get Schrenk's chord
c_schrenk_z = (c_z + c_elliptic_z) / 2;

I = sum(c_schrenk_z) * dz;
K = W / 2*I;


%% SPAR REQUIRED 1) IGNORING 2) CONSIDERING - SKIN MOMENT OF INERTIA



skin_thickness = 0.5 * 1e-3;
spar_thickness = 0.5 * 1e-3;

length_scalingfactor_z = Length_ScalingFactor_at_z(z, b, lambda);

I_xx = zeros(1, length(z));
I_yy = zeros(1, length(z));
I_xy = zeros(1, length(z));

% for i_z = 1:length(z)
%     x = x_root * length_scalingfactor_z;
%     y = y_root * length_scalingfactor_z;
%     
%     % Moments of Inertia of Skin at z
% %     [I_xx(i_z), I_yy(i_z), I_xy(i_z)] = Moments_of_Inertia(x, y, skin_thickness);
% 
% end



%% LOCATE SHEAR CENTER - SHEAR FLOW CALCULATIONS - 3 BOOM IDEALIZATION - CONSIDERING 1 SPAR AT MAXIMUM THICKNESS

% Set material properties
Aluminium_skin = Aluminium_7075();
Aluminium_spar = Aluminium_7075();

% Get the pre-requisite values required for shear flow calculations
% Discretization of the skin and spar
% ASSUMPTION: The spar endpoints are coincident with points on the skin
x = x_all_original;
y = y_all_original;
x_spar = x_spar_original;
y_spar = y_spar_original;


% Skin
panel_lengths = sqrt( (x([2:end, 1]) - x).^2 + (y([2:end, 1]) - y).^2 );
skin_thickness;
panel_areas = panel_lengths * skin_thickness;
midpoints_x = (x + x([2: end, 1])) / 2;
midpoints_y = (y + y([2: end, 1])) / 2;
% Spar
spar_panel_lengths = sqrt( (x_spar(2: end) - x_spar(1: end - 1)).^2 + (y_spar(2: end) - y_spar(1: end - 1)).^2 );
spar_thickness;
spar_panel_areas = spar_panel_lengths * spar_thickness;
spar_midpoints_x = (x_spar(1: end - 1) + x_spar(2: end)) / 2;
spar_midpoints_y = (y_spar(1: end - 1) + y_spar(2: end)) / 2;

% Find Centroid of section
centroid_x = (sum(midpoints_x .* panel_areas) + sum(spar_midpoints_x .* spar_panel_areas)) / (sum(panel_areas) + sum(spar_panel_areas));
centroid_y = (sum(midpoints_y .* panel_areas) + sum(spar_midpoints_y .* spar_panel_areas)) / (sum(panel_areas) + sum(spar_panel_areas));

% Get the coordinates of the points in the coordinate system with origin at centroid
x = x - centroid_x;
y = y - centroid_y;
x_spar = x_spar - centroid_x;
y_spar = y_spar - centroid_y;

% Indexes of points where the spar is located
spar_Idx_1;
spar_Idx_2;

[shearcenter_x, shearcenter_y] = ShearCenter_Location_func(x, y, skin_thickness, x_spar, y_spar, spar_thickness, spar_Idx_1, spar_Idx_2, Aluminium_skin, Aluminium_spar);

% Get Location of the shear center at the wing root
shearcenter_x_root = shearcenter_x * (c_root / chord_original);
shearcenter_y_root = shearcenter_y * (c_root / chord_original);

% Get other coordinates at the root
x_root = x * (c_root / chord_original);
y_root = y * (c_root / chord_original);
x_spar_root = x_spar * (c_root / chord_original);
y_spar_root = y_spar * (c_root / chord_original);

% figure(31);
% hold on; axis equal; grid on;
% plot(x_root, y_root);
% plot(x_spar_root, x_spar_root);
% plot(shearcenter_x_root, shearcenter_y_root, 'o');


%% SHEAR FLOW CALCULATIONS ALONG THE SPAN

% Get shear flow at z
z;
z = 0;

% Length Scaling factor w.r.t root
length_scalingfactor_z = Length_ScalingFactor_at_z(z, b, lambda);

% Create vectors for Moments of Inertia
I_xx = zeros(1, length(z));
I_yy = zeros(1, length(z));
I_xy = zeros(1, length(z));

V_x = zeros(1, length(z));
V_y = zeros(1, length(z));

for i_z = 1:length(z)

    % Coordinates of points on the skin and spar at z
    x = x_root * length_scalingfactor_z(i_z);
    y = y_root * length_scalingfactor_z(i_z);
    x_spar = x_spar_root * length_scalingfactor_z(i_z);
    y_spar = y_spar_root * length_scalingfactor_z(i_z);
    % Get location of shear center
    shearcenter_x = shearcenter_x_root * length_scalingfactor_z(i_z);
    shearcenter_y = shearcenter_y_root * length_scalingfactor_z(i_z);
    
    % Indexes of points where the spar is located
    spar_Idx_1;
    spar_Idx_2;

    % Moments of Area
    [I_xx(i_z), I_yy(i_z), I_xy(i_z)] = Moments_of_Inertia(x, y, skin_thickness, x_spar, y_spar, spar_thickness);
    D = I_xx(i_z) * I_yy(i_z) - I_xy(i_z)^2;
    
    % Location of the Shear Forces is assumed to be at quarter chord
%     V_location_x = c_quarter_x_root * length_scalingfactor_z(i_z);
%     V_location_y = c_quarter_y_root * length_scalingfactor_z(i_z);
    V_location_x = shearcenter_x;
    V_location_y = shearcenter_y;
    % Shear Forces and Torque at z
    V_x(i_z) = 0;
    V_y(i_z) = 51.3554;
    M_z(i_z) = 0;
    
    % Boom coordinates for skin
    N_boom = 2 * length(x);
    % NOTE: Number of panels after boom creation = N_boom
    % Odd points are boom split corners, even points are boom split centers
    boom_x = zeros(N_boom, 1);
    boom_x(1: 2: N_boom) = x;
    boom_x(2: 2: N_boom) = (x + x([2:end, 1])) / 2;
    boom_y = zeros(N_boom, 1);
    boom_y(1: 2: N_boom) = y;
    boom_y(2: 2: N_boom) = (y + y([2:end, 1])) / 2;
    % Boom coordinates for spar
    N_boom_spar = 2 * length(x_spar) - 1;
    spar_boom_x = zeros(N_boom_spar, 1);
    spar_boom_x(1: 2: N_boom_spar) = x_spar;
    spar_boom_x(2: 2: N_boom_spar) = (x_spar(1: end - 1) + x_spar(2: end)) / 2; 
    spar_boom_y = zeros(N_boom_spar, 1);
    spar_boom_y(1: 2: N_boom_spar) = y_spar;
    spar_boom_y(2: 2: N_boom_spar) = (y_spar(1: end - 1) + y_spar(2: end)) / 2; 
    % Get the boom indexes of points where the spar is located - these indices are along the airfoil
    spar_boom_Idx_1 = 2 * spar_Idx_1 - 1;
    spar_boom_Idx_2 = 2 * spar_Idx_2 - 1;

    % Get areas of the booms
    % Panel length = the length of current panel size. This is not the panel size after boom idealization
    panel_lengths = sqrt( (x([2:end, 1]) - x).^2 + (y([2:end, 1]) - y).^2 );
    skin_thickness;
    boom_areas = zeros(N_boom, 1);
    boom_areas(1: 2: N_boom) = (1/6) * panel_lengths * skin_thickness + (1/6) * panel_lengths([end, 1:end-1]) * skin_thickness;
    boom_areas(2: 2: N_boom) = (2/3) * panel_lengths * skin_thickness;
    % Boom areas for spar
    spar_panel_lengths = sqrt( (x_spar(2: end) - x_spar(1: end - 1)).^2 + (y_spar(2: end) - y_spar(1: end - 1)).^2 );
    spar_thickness;
    spar_boom_areas = zeros(N_boom_spar, 1);
    spar_boom_areas(1: 2: N_boom_spar) = [(1/6) * spar_panel_lengths * spar_thickness; 0] + [0; (1/6) * spar_panel_lengths * spar_thickness];
    spar_boom_areas(2: 2: N_boom_spar) = (2/3) * spar_panel_lengths * spar_thickness;
    % Add the spar boom areas at the 2 points where it joins with the skin
    boom_areas(spar_boom_Idx_1) = boom_areas(spar_boom_Idx_1) + spar_boom_areas(1);
    boom_areas(spar_boom_Idx_2) = boom_areas(spar_boom_Idx_2) + spar_boom_areas(end);
    % Change spar boom areas to the above values as well - for good measure
    spar_boom_areas(1) = boom_areas(spar_boom_Idx_1);
    spar_boom_areas(end) = boom_areas(spar_boom_Idx_2);

    % Get Matrix for each loop. Full skin loop, and, fore and aft loop.
    % Each matrix has columns as - boom_x, boom_y, boom area, shear flow(q), panel material shear modulus, panel thickness for that respective loop.
    Loop_full = [boom_x, boom_y, boom_areas, zeros(N_boom, 1)];
    % Create fore and aft loops with x, y, boom areas
    Loop_aft = [[boom_x(1: spar_boom_Idx_1); spar_boom_x(2: end-1); boom_x(spar_boom_Idx_2: end)], [boom_y(1: spar_boom_Idx_1); spar_boom_y(2: end-1); boom_y(spar_boom_Idx_2: end)], [boom_areas(1: spar_boom_Idx_1); spar_boom_areas(2: end-1); boom_areas(spar_boom_Idx_2: end)]];
    Loop_fore = [[boom_x(spar_boom_Idx_1: spar_boom_Idx_2); flip(spar_boom_x(2: end-1))], [boom_y(spar_boom_Idx_1: spar_boom_Idx_2); flip(spar_boom_y(2: end-1))], [boom_areas(spar_boom_Idx_1: spar_boom_Idx_2); flip(spar_boom_areas(2: end-1))]];
    
    % Determine loop lengths and get indices of critical points
    Loop_aft_length = length(Loop_aft(:, 1));
    Loop_aft_Idx_1 = spar_boom_Idx_1;
    Loop_aft_Idx_2 = spar_boom_Idx_1 + N_boom_spar - 1;
    Loop_fore_length = length(Loop_fore(:, 1));
    Loop_fore_Idx_1 = spar_boom_Idx_2 - spar_boom_Idx_1 + 1;

    % Add shear flow column to both
    Loop_aft = [Loop_aft, zeros(Loop_aft_length, 1)];
    Loop_fore = [Loop_fore, zeros(Loop_fore_length, 1)];
    % Add shear modulus of each panel for twist calculation
    Loop_aft = [Loop_aft, [ones(Loop_aft_Idx_1 - 1, 1) * Aluminium_skin.ShearModulus; ones(Loop_aft_Idx_2 - Loop_aft_Idx_1, 1) * Aluminium_spar.ShearModulus; ones(Loop_aft_length - Loop_aft_Idx_2 + 1, 1) * Aluminium_skin.ShearModulus]];
    Loop_fore = [Loop_fore, [ones(Loop_fore_Idx_1 - 1, 1) * Aluminium_skin.ShearModulus; ones(Loop_fore_length - Loop_fore_Idx_1 + 1, 1) * Aluminium_spar.ShearModulus]];
    % Add thickness of each panel for twist calculation
    Loop_aft = [Loop_aft, [ones(Loop_aft_Idx_1 - 1, 1) * skin_thickness; ones(Loop_aft_Idx_2 - Loop_aft_Idx_1, 1) * spar_thickness; ones(Loop_aft_length - Loop_aft_Idx_2 + 1, 1) * skin_thickness]];
    Loop_fore = [Loop_fore, [ones(Loop_fore_Idx_1 - 1, 1) * skin_thickness; ones(Loop_fore_length - Loop_fore_Idx_1 + 1, 1) * spar_thickness]];
    
    % Store 0 shear flow values at appropriate locations
    q0 = 0;
    q1 = 0;
    Loop_full(1, 4) = q0;
    Loop_aft(1, 4) = q0;
    Loop_aft(Loop_aft_Idx_1, 4) = q1;

    % Integrate in each loop to get the distribution of shear flow as a function of the unknown variables
    % Get the coefficients of shear flow formula
    coeff_1 = -(V_x(i_z) * I_xx(i_z) - V_y(i_z) * I_xy(i_z)) / D;
    coeff_2 = -(V_y(i_z) * I_yy(i_z) - V_x(i_z) * I_xy(i_z)) / D;

    % Integrate aft loop upto where spar joins skin
    Integration_Idx = 2: Loop_aft_Idx_1 - 1;
    dq = coeff_1 * Loop_aft(Integration_Idx, 1) .* Loop_aft(Integration_Idx, 3) + coeff_2 * Loop_aft(Integration_Idx, 2) .* Loop_aft(Integration_Idx, 3);
    Loop_aft(Integration_Idx, 4) = Loop_aft(Integration_Idx(1) - 1, 4) + cumsum(dq);
    % Set first panel of spar to q1
    Loop_aft(Loop_aft_Idx_1, 4) = q1;
    % Get shear flow along spar in aft loop
    Integration_Idx = Loop_aft_Idx_1 + 1: Loop_aft_Idx_2 - 1;
    dq = coeff_1 * Loop_aft(Integration_Idx, 1) .* Loop_aft(Integration_Idx, 3) + coeff_2 * Loop_aft(Integration_Idx, 2) .* Loop_aft(Integration_Idx, 3);
    Loop_aft(Integration_Idx, 4) = Loop_aft(Integration_Idx(1) - 1, 4) + cumsum(dq);
    
    % Integrate along fore loop
    % Set first panel
    Integration_Idx = 1;
    dq = coeff_1 * Loop_aft(Loop_aft_Idx_1, 1) .* Loop_aft(Loop_aft_Idx_1, 3) + coeff_2 * Loop_aft(Loop_aft_Idx_1, 2) .* Loop_aft(Loop_aft_Idx_1, 3);
    Loop_fore(Integration_Idx, 4) = Loop_aft(Loop_aft_Idx_1 - 1, 4) + dq - q1;
    % Get shear flow for rest of the fore loop
    Integration_Idx = 2: Loop_fore_Idx_1 - 1;
    dq = coeff_1 * Loop_fore(Integration_Idx, 1) .* Loop_fore(Integration_Idx, 3) + coeff_2 * Loop_fore(Integration_Idx, 2) .* Loop_fore(Integration_Idx, 3);
    Loop_fore(Integration_Idx, 4) = Loop_fore(Integration_Idx(1) - 1, 4) + cumsum(dq);
    % Set spar panels for fore loop = negative of those in aft loop
    Integration_Idx = Loop_fore_Idx_1: Loop_fore_length;
    Loop_fore(Integration_Idx, 4) = - flip(Loop_aft(Loop_aft_Idx_1: Loop_aft_Idx_2 - 1, 4));

    % Integrate along aft loop again
    % Get first panel
    Integration_Idx = Loop_aft_Idx_2;
    dq = coeff_1 * Loop_aft(Integration_Idx, 1) .* Loop_aft(Integration_Idx, 3) + coeff_2 * Loop_aft(Integration_Idx, 2) .* Loop_aft(Integration_Idx, 3);
    Loop_aft(Integration_Idx, 4) = Loop_aft(Integration_Idx - 1, 4) + Loop_fore(Loop_fore_Idx_1 - 1, 4) + dq;
    % Integrate along aft loop from spar top to trailing edge
    Integration_Idx = Loop_aft_Idx_2 + 1: Loop_aft_length;
    dq = coeff_1 * Loop_aft(Integration_Idx, 1) .* Loop_aft(Integration_Idx, 3) + coeff_2 * Loop_aft(Integration_Idx, 2) .* Loop_aft(Integration_Idx, 3);
    Loop_aft(Integration_Idx, 4) = Loop_aft(Integration_Idx(1) - 1, 4) + cumsum(dq);
    
    % Define variables for the unknown shear flows - define their location also
    syms q0 q1
    % Find the twist of each section
    % alpha = twist = integral ( q / (2 * mu * thickness * Area));
    A_aft = Area_enclosed_by_loop(Loop_aft(:, 1), Loop_aft(:, 2), Loop_aft(:, 6));
    A_fore = Area_enclosed_by_loop(Loop_fore(:, 1), Loop_fore(:, 2), Loop_fore(:, 6));
    
    % Get aft twist rate
    twist_aft = sum( Loop_aft(:, 4) ./ (2 * Loop_aft(:, 5) .* Loop_aft(:, 6) * A_aft) );
    Integration_Idx_aft_1 = 1: Loop_aft_Idx_1 - 1;
    Integration_Idx_aft_2 = Loop_aft_Idx_1: Loop_aft_Idx_2 - 1;
    Integration_Idx_aft_3 = Loop_aft_Idx_2: Loop_aft_length;
    twist_aft = twist_aft + sum( q0 ./ (2 * Loop_aft(Integration_Idx_aft_1, 5) .* Loop_aft(Integration_Idx_aft_1, 6) * A_aft) ) + sum( q1 ./ (2 * Loop_aft(Integration_Idx_aft_2, 5) .* Loop_aft(Integration_Idx_aft_2, 6) * A_aft) ) + sum( q0 ./ (2 * Loop_aft(Integration_Idx_aft_3, 5) .* Loop_aft(Integration_Idx_aft_3, 6) * A_aft) );
    % Get fore twist rate
    twist_fore = sum( Loop_fore(:, 4) ./ (2 * Loop_fore(:, 5) .* Loop_fore(:, 6) * A_fore) );
    Integration_Idx_fore_1 = 1: Loop_fore_Idx_1 - 1;
    Integration_Idx_fore_2 = Loop_fore_Idx_1: Loop_fore_length;
    twist_fore = twist_fore + sum( (q0 - q1) ./ (2 * Loop_fore(Integration_Idx_fore_1, 5) .* Loop_fore(Integration_Idx_fore_1, 6) * A_fore) ) + sum( (-q1) ./ (2 * Loop_fore(Integration_Idx_fore_2, 5) .* Loop_fore(Integration_Idx_fore_2, 6) * A_fore) );
    
    
    % To find shear center we assume that the shear forces pass through shear center. Therefore, twist of each section is 0.
    % Equate twist to 0 to get the unknown shear flows
    eqns = [twist_aft == 0, twist_fore == 0];
    Q = solve(eqns, [q0, q1]);
    q0 = double(Q.q0);
    q1 = double(Q.q1);
    
    % Get final shear flow distribution
    Loop_aft(Integration_Idx_aft_1, 4) = Loop_aft(Integration_Idx_aft_1, 4) + q0;
    Loop_aft(Integration_Idx_aft_2, 4) = Loop_aft(Integration_Idx_aft_2, 4) + q1;
    Loop_aft(Integration_Idx_aft_3, 4) = Loop_aft(Integration_Idx_aft_3, 4) + q0;
    Loop_fore(Integration_Idx_fore_1, 4) = Loop_fore(Integration_Idx_fore_1, 4) + q0 - q1;
    Loop_fore(Integration_Idx_fore_2, 4) = Loop_fore(Integration_Idx_fore_2, 4) + (-q1);

    % Consider torque due to shear forces about the shear center
    M_z_Shearforces = - V_y(i_z) * (shearcenter_x - V_location_x) + V_x(i_z) * (shearcenter_y - V_location_y);
    M_z_Moment = M_z(i_z);
    M_z_total = M_z_Shearforces + M_z_Moment;
    
    % Assume some unknown constant shear flow in each loop
    syms q_aft_torque q_fore_torque
    A_aft;
    A_fore;
    % Equate Torque applied to torque by assumed shear flow distribution
    % NOTE: Consider sign convention
    eqn1 = M_z_total == - (q_aft_torque * A_aft + q_fore_torque * A_fore);
    
    % Find the twist angles due to the assumed shear flow distributions
    Idx_1 = 1: (Loop_aft_Idx_1 - 1);
    Idx_2 = Loop_aft_Idx_1: (Loop_aft_Idx_2 - 1);
    Idx_3 = Loop_aft_Idx_2: Loop_aft_length;
    twist_aft = sum( q_aft_torque ./ (2 * Loop_aft([Idx_1, Idx_3], 5) .* Loop_aft([Idx_1, Idx_3], 6) * A_aft) ) + sum( (q_aft_torque - q_fore_torque) ./ (2 * Loop_aft(Idx_2, 5) .* Loop_aft(Idx_2, 6) * A_aft) );
    Idx_1 = 1: (Loop_fore_Idx_1 - 1);
    Idx_2 = Loop_fore_Idx_1: Loop_fore_length;
    twist_fore = sum( q_fore_torque ./ (2 * Loop_fore(Idx_1, 5) .* Loop_fore(Idx_1, 6) * A_fore) ) + sum( (q_fore_torque - q_aft_torque) ./ (2 * Loop_fore(Idx_2, 5) .* Loop_fore(Idx_2, 6) * A_fore) );
    
    % Compatibility condition -> twist angles are equal
    eqn2 = twist_aft == twist_fore;
    
    % Solve to get shear flow due to torque
    eqns = [eqn1, eqn2];
    Q = solve(eqns, [q_aft_torque, q_fore_torque]);
    q_aft_torque = double(Q.q_aft_torque);
    q_fore_torque = double(Q.q_fore_torque);
    
    % Superimpose shear flow due to shear force and torque to get final shear flow distribution
    Loop_aft(:, 4) = Loop_aft(:, 4) + q_aft_torque;
    Loop_aft((Loop_aft_Idx_1: (Loop_aft_Idx_2 - 1)), 4) = Loop_aft((Loop_aft_Idx_1: (Loop_aft_Idx_2 - 1)), 4) - q_fore_torque;
    Loop_fore(:, 4) = Loop_fore(:, 4) + q_fore_torque;
    Loop_fore((Loop_fore_Idx_1: Loop_fore_length), 4) = Loop_fore((Loop_fore_Idx_1: Loop_fore_length), 4) - q_aft_torque;


end



figure(10);
hold on; grid on;
a = 1: N_boom;
y = [Loop_aft(1:(Loop_aft_Idx_1-1), 4); Loop_fore(1:(Loop_fore_Idx_1-1), 4); Loop_aft(Loop_aft_Idx_2: end, 4)];
plot(a, y);

% return
% plot(1:length(Loop_aft(1:(Loop_aft_Idx_1-1), 1)), Loop_aft(1:(Loop_aft_Idx_1-1), 4));
% plot(1:length(Loop_fore(1:(Loop_fore_Idx_1-1), 1)), Loop_fore(1:(Loop_fore_Idx_1-1), 4));
% plot(1:length(Loop_aft(Loop_aft_Idx_2: end, 1)), Loop_aft(Loop_aft_Idx_2: end, 4));






























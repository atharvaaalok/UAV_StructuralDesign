function [ShearCenter_x, ShearCenter_y] = ShearCenter_Location_func(x, y, skin_thickness, x_spar, y_spar, spar_thickness, spar_Idx_1, spar_Idx_2, Aluminium_skin, Aluminium_spar)

    % Moments of Area
    [I_xx, I_yy, I_xy] = Moments_of_Inertia(x, y, skin_thickness, x_spar, y_spar, spar_thickness);

    D = I_xx * I_yy - I_xy^2;
    

    %% FIND SHEAR CENTER X COORDINATE

    % Appy only V_y for finding the x-coordinate of the shear center
    V_x = 0;
    V_y = .1;
    
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
    
%     figure(3);
%     hold on; grid on; axis equal
%     plot(x, y);
%     plot(boom_x, boom_y, 'o');
%     plot(spar_boom_x, spar_boom_y, 'o');
%     plot(boom_x(spar_boom_Idx_1), boom_y(spar_boom_Idx_1), 'Marker', 'x');
%     plot(boom_x(spar_boom_Idx_2), boom_y(spar_boom_Idx_2), 'Marker', 'x');

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
%     % Plot the points to confirm
%     figure(4);
%     hold on; grid on; axis equal
%     plot(x, y);
%     plot(boom_x, boom_y, 'o');
%     plot(spar_boom_x, spar_boom_y, 'o');
%     plot(Loop_aft(Loop_aft_Idx_1, 1), Loop_aft(Loop_aft_Idx_1, 2), 'Marker', 'x');
%     plot(Loop_aft(Loop_aft_Idx_2, 1), Loop_aft(Loop_aft_Idx_2, 2), 'Marker', 'x');
%     figure(5);
%     hold on; grid on; axis equal
%     plot(x, y);
%     plot(boom_x, boom_y, 'o');
%     plot(spar_boom_x, spar_boom_y, 'o');
%     plot(Loop_fore(1, 1), Loop_fore(1, 2), 'Marker', 'x');
%     plot(Loop_fore(Loop_fore_Idx_1, 1), Loop_fore(Loop_fore_Idx_1, 2), 'Marker', 'x');

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
    coeff_1 = -(V_x * I_xx - V_y * I_xy) / D;
    coeff_2 = -(V_y * I_yy - V_x * I_xy) / D;

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
    
    % Define unknown variables for shear center coordinates
    % Get moment of Shear forces at the centroid assume they act at shear center
    % Follow sign convention
    syms shearcenter_x
    M_z_Shearforces = V_y * shearcenter_x;
    
    % Get shear flow direction vector for aft loop
    shearflow_direction_aft = [(Loop_aft([2: end, 1], 1) - Loop_aft(:, 1)), (Loop_aft([2: end, 1], 2) - Loop_aft(:, 2))];
    shearflow_direction_aft = shearflow_direction_aft ./ (sqrt(shearflow_direction_aft(:, 1).^2 + shearflow_direction_aft(:, 2).^2));
    shearflow_direction_aft(isnan(shearflow_direction_aft)) = 0;
    shearflow_vector_aft = Loop_aft(:, 4) .*  shearflow_direction_aft;
    % Find shear flow vector for fore loop without spar
    shearflow_direction_fore = [(Loop_fore(2: Loop_fore_Idx_1, 1) - Loop_fore(1: Loop_fore_Idx_1 - 1, 1)), (Loop_fore(2: Loop_fore_Idx_1, 2) - Loop_fore(1: Loop_fore_Idx_1 - 1, 2))];
    shearflow_direction_fore = shearflow_direction_fore ./ (sqrt(shearflow_direction_fore(:, 1).^2 + shearflow_direction_fore(:, 2).^2));
    shearflow_direction_fore(isnan(shearflow_direction_fore)) = 0;
    shearflow_vector_fore = Loop_fore(1: Loop_fore_Idx_1 - 1, 4) .*  shearflow_direction_fore;
    % Get panel lengths
    panel_lengths_aft = sqrt( (Loop_aft([2:end, 1], 1) - Loop_aft(:, 1)).^2 + (Loop_aft([2:end, 1], 2) - Loop_aft(:, 2)).^2 );
    panel_lengths_fore = sqrt( (Loop_fore(2: Loop_fore_Idx_1, 1) - Loop_fore(1: Loop_fore_Idx_1 - 1, 1)).^2 + (Loop_fore(2: Loop_fore_Idx_1, 2) - Loop_fore(1: Loop_fore_Idx_1 - 1, 2)).^2 );
    
    % Get Shear Force vectors
    shearforce_vector_aft = shearflow_vector_aft .* panel_lengths_aft;
    shearforce_vector_fore = shearflow_vector_fore .* panel_lengths_fore;
    
    % Cross product Moment = r x F
    % Get moment of the shear flow at the centroid
    % Get direction vector of panel midpoint from centroid
    midpoint_vector_aft = [(Loop_aft([2: end, 1], 1) + Loop_aft(:, 1)) / 2, (Loop_aft([2: end, 1], 2) + Loop_aft(:, 2)) / 2];
    midpoint_vector_fore = [(Loop_fore(1: Loop_fore_Idx_1 - 1, 1) + Loop_fore(2: Loop_fore_Idx_1, 1)) / 2, (Loop_fore(1: Loop_fore_Idx_1 - 1, 2) + Loop_fore(2: Loop_fore_Idx_1, 2)) / 2];
    % For cross product 3 elements needed in each row
    shearforce_vector_aft = [shearforce_vector_aft, zeros(size(shearforce_vector_aft, 1), 1)];
    shearforce_vector_fore = [shearforce_vector_fore, zeros(size(shearforce_vector_fore, 1), 1)];
    midpoint_vector_aft = [midpoint_vector_aft, zeros(size(midpoint_vector_aft, 1), 1)];
    midpoint_vector_fore = [midpoint_vector_fore, zeros(size(midpoint_vector_fore, 1), 1)];

    % Find cross product
    M_z_shearflow_aft = cross(midpoint_vector_aft, shearforce_vector_aft);
    M_z_shearflow_aft = M_z_shearflow_aft(:, 3);
    M_z_shearflow_fore = cross(midpoint_vector_fore, shearforce_vector_fore);
    M_z_shearflow_fore = M_z_shearflow_fore(:, 3);
    % Get total moment due to shear flow
    M_z_shearflow = sum(M_z_shearflow_aft) + sum(M_z_shearflow_fore);
    M_z_shearflow = double(M_z_shearflow);
    
    % Equate to get the shear center
    eqn = M_z_Shearforces == M_z_shearflow;
    shearcenter_x = solve(eqn, shearcenter_x);
    shearcenter_x = double(shearcenter_x);
    

    %% FIND SHEAR CENTER Y COORDINATE

    % Shear Forces and Torque
    V_x = 1;
    V_y = 0;
    
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
    
%     figure(3);
%     hold on; grid on; axis equal
%     plot(x, y);
%     plot(boom_x, boom_y, 'o');
%     plot(spar_boom_x, spar_boom_y, 'o');
%     plot(boom_x(spar_boom_Idx_1), boom_y(spar_boom_Idx_1), 'Marker', 'x');
%     plot(boom_x(spar_boom_Idx_2), boom_y(spar_boom_Idx_2), 'Marker', 'x');

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
    % Plot the points to confirm
    figure(4);
    hold on; grid on; axis equal
    plot(x, y);
    plot(boom_x, boom_y, 'o');
    plot(spar_boom_x, spar_boom_y, 'o');
    plot(Loop_aft(Loop_aft_Idx_1, 1), Loop_aft(Loop_aft_Idx_1, 2), 'Marker', 'x');
    plot(Loop_aft(Loop_aft_Idx_2, 1), Loop_aft(Loop_aft_Idx_2, 2), 'Marker', 'x');
    figure(5);
    hold on; grid on; axis equal
    plot(x, y);
    plot(boom_x, boom_y, 'o');
    plot(spar_boom_x, spar_boom_y, 'o');
    plot(Loop_fore(1, 1), Loop_fore(1, 2), 'Marker', 'x');
    plot(Loop_fore(Loop_fore_Idx_1, 1), Loop_fore(Loop_fore_Idx_1, 2), 'Marker', 'x');

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
    coeff_1 = -(V_x * I_xx - V_y * I_xy) / D;
    coeff_2 = -(V_y * I_yy - V_x * I_xy) / D;
    
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

    % Define unknown variables for shear center coordinates
    % Get moment of Shear forces at the centroid assume they act at shear center
    % Follow sign convention
    syms shearcenter_y
    M_z_Shearforces = - V_x * shearcenter_y;
    
    % Get shear flow direction vector for aft loop
    shearflow_direction_aft = [(Loop_aft([2: end, 1], 1) - Loop_aft(:, 1)), (Loop_aft([2: end, 1], 2) - Loop_aft(:, 2))];
    shearflow_direction_aft = shearflow_direction_aft ./ (sqrt(shearflow_direction_aft(:, 1).^2 + shearflow_direction_aft(:, 2).^2));
    shearflow_direction_aft(isnan(shearflow_direction_aft)) = 0;
    shearflow_vector_aft = Loop_aft(:, 4) .*  shearflow_direction_aft;
    % Find shear flow vector for fore loop without spar
    shearflow_direction_fore = [(Loop_fore(2: Loop_fore_Idx_1, 1) - Loop_fore(1: Loop_fore_Idx_1 - 1, 1)), (Loop_fore(2: Loop_fore_Idx_1, 2) - Loop_fore(1: Loop_fore_Idx_1 - 1, 2))];
    shearflow_direction_fore = shearflow_direction_fore ./ (sqrt(shearflow_direction_fore(:, 1).^2 + shearflow_direction_fore(:, 2).^2));
    shearflow_direction_fore(isnan(shearflow_direction_fore)) = 0;
    shearflow_vector_fore = Loop_fore(1: Loop_fore_Idx_1 - 1, 4) .*  shearflow_direction_fore;
    % Get panel lengths
    panel_lengths_aft = sqrt( (Loop_aft([2:end, 1], 1) - Loop_aft(:, 1)).^2 + (Loop_aft([2:end, 1], 2) - Loop_aft(:, 2)).^2 );
    panel_lengths_fore = sqrt( (Loop_fore(2: Loop_fore_Idx_1, 1) - Loop_fore(1: Loop_fore_Idx_1 - 1, 1)).^2 + (Loop_fore(2: Loop_fore_Idx_1, 2) - Loop_fore(1: Loop_fore_Idx_1 - 1, 2)).^2 );
    
    % Get Shear Force vectors
    shearforce_vector_aft = shearflow_vector_aft .* panel_lengths_aft;
    shearforce_vector_fore = shearflow_vector_fore .* panel_lengths_fore;
    
    % Cross product Moment = r x F
    % Get moment of the shear flow at the centroid
    % Get direction vector of panel midpoint from centroid
    midpoint_vector_aft = [(Loop_aft([2: end, 1], 1) + Loop_aft(:, 1)) / 2, (Loop_aft([2: end, 1], 2) + Loop_aft(:, 2)) / 2];
    midpoint_vector_fore = [(Loop_fore(1: Loop_fore_Idx_1 - 1, 1) + Loop_fore(2: Loop_fore_Idx_1, 1)) / 2, (Loop_fore(1: Loop_fore_Idx_1 - 1, 2) + Loop_fore(2: Loop_fore_Idx_1, 2)) / 2];
    % For cross product 3 elements needed in each row
    shearforce_vector_aft = [shearforce_vector_aft, zeros(size(shearforce_vector_aft, 1), 1)];
    shearforce_vector_fore = [shearforce_vector_fore, zeros(size(shearforce_vector_fore, 1), 1)];
    midpoint_vector_aft = [midpoint_vector_aft, zeros(size(midpoint_vector_aft, 1), 1)];
    midpoint_vector_fore = [midpoint_vector_fore, zeros(size(midpoint_vector_fore, 1), 1)];

    % Find cross product
    M_z_shearflow_aft = cross(midpoint_vector_aft, shearforce_vector_aft);
    M_z_shearflow_aft = M_z_shearflow_aft(:, 3);
    M_z_shearflow_fore = cross(midpoint_vector_fore, shearforce_vector_fore);
    M_z_shearflow_fore = M_z_shearflow_fore(:, 3);
    % Get total moment due to shear flow
    M_z_shearflow = sum(M_z_shearflow_aft) + sum(M_z_shearflow_fore);
    M_z_shearflow = double(M_z_shearflow);
    
    % Equate to get the shear center
    eqn = M_z_Shearforces == M_z_shearflow;
    shearcenter_y = solve(eqn, shearcenter_y);
    shearcenter_y = double(shearcenter_y);


    %% SET FINAL VALUES OF THE COORDINATES

    ShearCenter_x = shearcenter_x;
    ShearCenter_y = shearcenter_y;
    
    PS = PLOT_STANDARDS();

    fig1_comps.fig = figure(31);
    hold on; axis equal; grid on;
    fig1_comps.p1 = plot(x, y, 'LineWidth', 1.75);
    fig1_comps.p2 = plot(x_spar, y_spar, 'LineWidth', 3);
    fig1_comps.p3 = plot(0, 0, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', PS.Green1, 'MarkerEdgeColor', PS.Green3);
    fig1_comps.p4 = plot(ShearCenter_x, ShearCenter_y, 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 10, 'MarkerFaceColor', PS.Orange1, 'MarkerEdgeColor', PS.Orange2);

    STANDARDIZE_FIGURE(fig1_comps);
    SAVE_MY_FIGURE(fig1_comps, 'Figures/ShearCenter_CG_Locations.png', 'small');



    double(ShearCenter_x)
    double(ShearCenter_y)




end
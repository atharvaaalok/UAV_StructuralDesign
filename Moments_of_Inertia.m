function [I_xx, I_yy, I_xy] = Moments_of_Inertia(x, y, skin_thickness, x_spar, y_spar, spar_thickness)

    % Skin
    panel_lengths = sqrt( (x([2:end, 1]) - x).^2 + (y([2:end, 1]) - y).^2 );
    panel_areas = panel_lengths * skin_thickness;
    midpoints_x = (x + x([2: end, 1])) / 2;
    midpoints_y = (y + y([2: end, 1])) / 2;
    % Spar
    spar_panel_lengths = sqrt( (x_spar(2: end) - x_spar(1: end - 1)).^2 + (y_spar(2: end) - y_spar(1: end - 1)).^2 );
    spar_panel_areas = spar_panel_lengths * spar_thickness;
    spar_midpoints_x = (x_spar(1: end - 1) + x_spar(2: end)) / 2;
    spar_midpoints_y = (y_spar(1: end - 1) + y_spar(2: end)) / 2;

    I_xx = sum(midpoints_y.^2 .* panel_areas) + sum(spar_midpoints_y.^2 .* spar_panel_areas);
    I_yy = sum(midpoints_x.^2 .* panel_areas) + sum(spar_midpoints_x.^2 .* spar_panel_areas);
    I_xy = sum(midpoints_x .* midpoints_y .* panel_areas) + sum(spar_midpoints_x .* spar_midpoints_y .* spar_panel_areas);


end






























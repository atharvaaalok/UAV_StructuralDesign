function Total_Area = Area_enclosed_by_loop(x, y, thickness_distribution)
    
    % Skin
    panel_lengths = sqrt( (x([2:end, 1]) - x).^2 + (y([2:end, 1]) - y).^2 );
    panel_areas = panel_lengths .* thickness_distribution;
    midpoints_x = (x + x([2: end, 1])) / 2;
    midpoints_y = (y + y([2: end, 1])) / 2;
    
    % Find Centroid of section
    centroid_x = sum(midpoints_x .* panel_areas) / sum(panel_areas);
    centroid_y = sum(midpoints_y .* panel_areas) / sum(panel_areas);
    
    x1 = x;
    x2 = x([2: end, 1]);
    x3 = x * 0 + centroid_x;
    y1 = y;
    y2 = y([2: end, 1]);
    y3 = y * 0 + centroid_y;

    Triangle_areas = abs( (1/2) * ( x1 .* (y2 - y3) + x2 .* (y3 - y1) + x3 .* (y1 - y2) ) );

    Total_Area = sum(Triangle_areas);
    
end
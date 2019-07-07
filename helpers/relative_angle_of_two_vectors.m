function [ angle_in_degree ] = relative_angle_of_two_vectors(vec1_x, vec1_y, vec2_x, vec2_y, vec1_icx, vec1_icy, vec2_icx, vec2_icy)
%func angle_of_two_vectors()
%   Angle is calculated in screen coordination. Return value is in degrees.
    ang1 = calculate_vector_angle(vec1_x, vec1_y, vec1_icx, vec1_icy);
    ang2 = calculate_vector_angle(vec2_x, vec2_y, vec2_icx, vec2_icy);
    ang_diff = ang1 - ang2;
    angle_in_degree = abs(ang_diff);
    
    if angle_in_degree > 180
        angle_in_degree = 360 - angle_in_degree;
    end
end


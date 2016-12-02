function [c] = correcting_surveys(Easting, Northing, OthoHgt)

if nargin < 3 % Run a demo
    Easting = [0.0882 0.0851 0.0822 0.0616 0.001 0.2132 0.254 0.4606 0.6791 0.9556];
    Northing = [27.1078 27.1108 27.1098 27.1078 27.0786 27.5766 28.1009 28.5976 29.0495 29.4665];
    OthoHgt = [0.333 0.3332 0.3338 0.3338 0.4059 0.3224 0.1945 0.0971 -0.018 -0.1295];
end

[w, h] = size(Easting);
if w > h
    Easting = Easting';
    Northing = Northing';
    OthoHgt = OthoHgt';
end

% Calculate the distance between the 1st & last northing & easting measurements
% Add 100 to both the raw easting and raw northing data so that errors are
% not induced by negative values
adjusted_easting = Easting + 100;
adjusted_northing = Northing + 100;

length_data = length(Easting);
extracted_easting (1,1) = adjusted_easting(1,1);
extracted_northing (1,1) = adjusted_northing(1,1);
extracted_easting(2,1) = adjusted_easting(1,length_data);
extracted_northing (2,1) = adjusted_northing(1,length_data);

part1 = diff(extracted_northing);
total_distance_northing = abs(part1);

part1 = diff(extracted_easting);
total_distance_easting = abs(part1);

% First, calculate the straight line distance between the 1st and last
% measurements. This is neccesary to compare the actual distance with the
% corrected distance at the end

% Calculate the angle a
part1 = total_distance_easting/total_distance_northing;
a_angle = atand(part1);
clear part1 

% Calculate the remaining angle b
total_profile_b_angle = 180 - a_angle - 90;

% Calculate the length of the cross section survey
part1 = sind (total_profile_b_angle);
cross_section_length = total_distance_northing/part1;
clear part1 a_angle b_angle length_data

% Assign the cells as zero due to the future points being relative to the
% first

a_angle (1,1) = 0;
b_angle (1,1) = 0;
x (1,1) = 0;
c_angle (1,1) = 0;
c (1,1) = 0;


% Calculate & change the relative points into absolute values
relative_easting = adjusted_easting - adjusted_easting(1,1);
relative_northing = adjusted_northing - adjusted_northing(1,1);

z = 2;
maxI = size (relative_easting) + 1;

while z < maxI (1,2)
    
    % Calculate the unknown angles of the outside triangle
    part1 = relative_northing (1,z) / relative_easting (1,z);
    a_angle (1,z) = atand(part1);
    a_angle_abs(1,z) = abs(a_angle(1,z));
    b_angle (1,z) = 180 - 90 - a_angle (1,z);

    
    % Calculate the distance between the first measurement and the current one
    part1 = sind(a_angle_abs(1,z));
    x (1,z) = relative_northing(1,z)/part1;


    % Calculate the remaining neccesary angles of the internal triangle
    c_angle (1,z) = 90 - b_angle (1,z);
    c_angle_abs (1,z) = abs(c_angle(1,z));
    d_angle (1,z) = total_profile_b_angle(1,1) - a_angle_abs(1,z);
    d_angle_abs (1,z) = abs(d_angle(1,z));
    x_angle (1,z) = 180 - d_angle_abs(1,z) - c_angle_abs(1,z);
    
    % Using the law of sines, calculate the length of the remaining sides
    part1 = x(1,z) .* sind(d_angle_abs(1,z)); 
    part2 = sind(x_angle(1,z));
    d(1,z) = part1 / part2;
    
    part1 = x(1,z) .* sind(c_angle_abs(1,z));
    part2 = sind(x_angle(1,z));
    c(1,z) = part1 / part2;
    c(1,z) = abs(c(1,z));
    
    % c (or ans in the demo) is the straight line distance for each survey point along transect
    % line i.e. survey points pulled into the transect line
    
    z = z + 1;

end

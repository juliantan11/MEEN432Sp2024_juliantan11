% Define vehicle parameters
Car_Velocity = 60;
sim_time = 1;

% Define track parameters
straight_length = 900;    % Length of straight sections in meters
curve_radius = 200;       % Radius of curved sections in meters
track_width = 15;         % Width of the track in meters
num_waypoints = 100;      % Number of waypoints

% Calculate delta_s for equal spacing of waypoints
total_length = 2 * straight_length + 2 * pi * curve_radius;
delta_s = total_length / (num_waypoints - 1);

% Calculate delta_theta for curved sections
delta_theta = delta_s / curve_radius;  

% Initialize arrays for coordinates
x = zeros(1, num_waypoints);
y = zeros(1, num_waypoints);
theta = zeros(1, num_waypoints);

% Generate waypoints for the track
for i = 1:num_waypoints
    s = (i - 1) * delta_s;
    
    % First straight section
    if s <= straight_length 
        x(i) = s;
        y(i) = 0;
        theta(i) = 0; % Straight section, no change in theta

    % First curved section
    elseif s <= straight_length + pi * curve_radius 
        segment_s = s - straight_length;
        theta(i) = segment_s / curve_radius;
        [x(i), y(i)] = rotate_point(straight_length, 0, straight_length, curve_radius, theta(i));
        
    % Second straight section
    elseif s <= 2 * straight_length + pi * curve_radius 
        segment_s = s - straight_length - pi * curve_radius;
        x(i) = straight_length - segment_s;
        y(i) = 2 * curve_radius;
        theta(i) = pi; % Straight section, no change in theta
        
    % Second curved section
    else 
        segment_s = s - 2 * straight_length - pi * curve_radius;
        theta(i) = pi + segment_s / curve_radius;
        [x(i), y(i)] = rotate_point(0, 0, 0, curve_radius, theta(i));
    end
end

% Call raceStat.m script
path.radius = 200; % Radius of Curves
path.width = 15;   % Width of the Track
path.l_st = 900;   % Length of Straightaways

Wp = [transpose(x), transpose(y)];
delta_f = 0;

% Function to rotate points around a center
function [x_rotated, y_rotated] = rotate_point(x, y, x_center, y_center, theta)
    R = [cos(theta), -sin(theta);
         sin(theta), cos(theta)];

    % Translate the point to be rotated to the origin
    x_translated = x - x_center;
    y_translated = y - y_center;

    % Perform the rotation using rotation matrix
    rotated_coords = R * [x_translated; y_translated];

    % Translate the rotated point back to its original position
    x_rotated = rotated_coords(1,:) + x_center;
    y_rotated = rotated_coords(2,:) + y_center;
end

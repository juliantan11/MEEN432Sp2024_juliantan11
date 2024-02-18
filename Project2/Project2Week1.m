% Track Specs
straightaway_length = 900; % Length of Track Straightaways [m]
turn_radius = 200; % Radius of Track Turns [m]
turn_angle = 180; % Angle of Turn (From Entry to Exit) [deg]
track_width = 15; % Width of Track [m]
car_width = 1; % Width of Car [m]
car_length = 2.5; % Length of Car (From Center) [m]

% Waypoint Setup
num_waypoints = 120; % Number of Waypoints for Track
delta_s = (2*straightaway_length)/(num_waypoints/2); % Change in Distance in Track Straightaways (30 for 120) [m]
delta_theta = (2*turn_angle)/(num_waypoints/2); % Change in Angle in Track Turns (6 for 120) [deg]
disp(delta_s)
disp(delta_theta)

% Track Coordinates
x_track_array = zeros(1,num_waypoints); % X Track Coordinates (Updated from 0 as code runs)
y_track_array = zeros(1,num_waypoints); % Y Track Coordinates (Updated from 0 as code runs)
theta_track_array = zeros(1,num_waypoints); % Angle Track Coordinates (Updated from 0 as code runs)

% Car Coordinates
car_RF_x_array = zeros(1,num_waypoints); % Car RF X Coordinates (Updated from 0 as code runs)
car_RR_x_array = zeros(1,num_waypoints); % Car RR X Coordinates (Updated from 0 as code runs)
car_LF_x_array = zeros(1,num_waypoints); % Car LF X Coordinates (Updated from 0 as code runs)
car_LR_x_array = zeros(1,num_waypoints); % Car LR X Coordinates (Updated from 0 as code runs)
car_RF_y_array = zeros(1,num_waypoints); % Car RF Y Coordinates (Updated from 0 as code runs)
car_RR_y_array = zeros(1,num_waypoints); % Car RR Y Coordinates (Updated from 0 as code runs)
car_LF_y_array = zeros(1,num_waypoints); % Car LF Y Coordinates (Updated from 0 as code runs)
car_LR_y_array = zeros(1,num_waypoints); % Car LR Y Coordinates (Updated from 0 as code runs)

for i = 1:num_waypoints
    % Front Straightaway
    if i <= num_waypoints/4
        x_track_array(i) = (i - 1)*delta_s;
        y_track_array(i) = 0;
        theta_track_array(i) = 0;
        car_RF_x_array(i) = (i - 1)*delta_s + car_length;
        car_RR_x_array(i) = (i - 1)*delta_s - car_length;
        car_LF_x_array(i) = (i - 1)*delta_s + car_length;
        car_LR_x_array(i) = (i - 1)*delta_s - car_length;
        car_RF_y_array(i) = -car_width;
        car_RR_y_array(i) = -car_width;
        car_LF_y_array(i) = car_width;
        car_LR_y_array(i) = car_width;
    % Turns 1 & 2
    elseif (num_waypoints/4 < i) && (i <= num_waypoints/2)
        theta_track_array(i) = (i - num_waypoints/4 - 1) * delta_theta;
        x_track_array(i) = x_track_array(num_waypoints/4) + turn_radius * sin(theta_track_array(i));
        y_track_array(i) = turn_radius * (1 - cos(theta_track_array(i)));
        car_RF_x_array(i) = x_track_array(num_waypoints/4) + (turn_radius + car_width) * sin(theta_track_array(i));
        car_RR_x_array(i) = x_track_array(num_waypoints/4) + (turn_radius + car_width) * sin(theta_track_array(i));
        car_LF_x_array(i) = x_track_array(num_waypoints/4) + (turn_radius - car_width) * sin(theta_track_array(i));
        car_LR_x_array(i) = x_track_array(num_waypoints/4) + (turn_radius - car_width) * sin(theta_track_array(i));
        car_RF_y_array(i) = (turn_radius + car_width) * (1 - cos(theta_track_array(i)));
        car_RR_y_array(i) = (turn_radius + car_width) * (1 - cos(theta_track_array(i)));
        car_LF_y_array(i) = (turn_radius - car_width) * (1 - cos(theta_track_array(i)));
        car_LR_y_array(i) = (turn_radius - car_width) * (1 - cos(theta_track_array(i)));
    % Back Straightaway
    elseif (num_waypoints/2 < i) && (i <= num_waypoints*(3/4))
        x_track_array(i) = (i - 1)*-delta_s;
        y_track_array(i) = 2*turn_radius;
        theta_track_array(i) = 0;
        car_RF_x_array(i) = (i - 1)*delta_s + car_length;
        car_RR_x_array(i) = (i - 1)*delta_s - car_length;
        car_LF_x_array(i) = (i - 1)*delta_s + car_length;
        car_LR_x_array(i) = (i - 1)*delta_s - car_length;
        car_RF_y_array(i) = 2*turn_radius + car_width;
        car_RR_y_array(i) = 2*turn_radius + car_width;
        car_LF_y_array(i) = 2*turn_radius - car_width;
        car_LR_y_array(i) = 2*turn_radius - car_width;
    % Turns 3 & 4
    else
        theta_track_array(i) = (i - num_waypoints * (3 / 4) - 1) * delta_theta;
        x_track_array(i) = x_track_array(num_waypoints * (3 / 4)) + turn_radius * sin(theta_track_array(i));
        y_track_array(i) = turn_radius - turn_radius * cos(theta_track_array(i));
        car_RF_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius + car_width) * sin(theta_track_array(i));
        car_RR_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius + car_width) * sin(theta_track_array(i));
        car_LF_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius - car_width) * sin(theta_track_array(i));
        car_LR_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius - car_width) * sin(theta_track_array(i));
        car_RF_y_array(i) = 2 * turn_radius - (turn_radius + car_width) * cos(theta_track_array(i));
        car_RR_y_array(i) = 2 * turn_radius - (turn_radius + car_width) * cos(theta_track_array(i));
        car_LF_y_array(i) = 2 * turn_radius - (turn_radius - car_width) * cos(theta_track_array(i));
        car_LR_y_array(i) = 2 * turn_radius - (turn_radius - car_width) * cos(theta_track_array(i));
    end
end

% Plotting the Track
figure;
plot(x_track_array, y_track_array, 'k', 'LineWidth', track_width);
hold on;

for i = 1:num_waypoints
    % Plotting the Car
    car_x_array = [car_RF_x_array(i),car_RR_x_array(i),car_LF_x_array(i),car_LR_x_array(i)];
    car_y_array = [car_RF_y_array(i),car_RR_y_array(i),car_LF_y_array(i),car_LR_y_array(i)];
    patch(car_x_array,car_y_array,'b')
    
    % Plotting the Car Path
    plot(x_track_array, y_track_array, 'r--', 'LineWidth', 1);
end

% Set axis limits
axis equal;
axis([-100, straightaway_length + 100, -50, 2 * turn_radius + 50]);

% Add labels and title
xlabel('X (m)');
ylabel('Y (m)');
title('Race Track with Racecar');

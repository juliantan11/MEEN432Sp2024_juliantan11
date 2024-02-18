% Track Specs
straightaway_length = 900; % Length of Track Straightaways [m]
turn_radius = 200; % Radius of Track Turns [m]
turn_angle = pi; % Angle of Turn (From Entry to Exit) [rad]
track_width = 15; % Width of Track [m]
car_width = 1; % Width of Car [m]
car_length = 2.5; % Length of Car (From Center) [m]

% Waypoint Setup
num_waypoints = 120; % Number of Waypoints for Track
delta_s = (2*straightaway_length)/(num_waypoints/2); % Change in Distance in Track Straightaways (30 for 120) [m]
delta_theta = (2*turn_angle)/(num_waypoints/2); % Change in Angle in Track Turns (pi/30 for 120) [rad]

% Track Coordinates
x_track_array = zeros(1,num_waypoints+1); % X Track Coordinates (Updated from 0 as code runs)
y_track_array = zeros(1,num_waypoints+1); % Y Track Coordinates (Updated from 0 as code runs)
theta_track_array = zeros(1,num_waypoints+1); % Angle Track Coordinates (Updated from 0 as code runs)

% Car Coordinates
car_RF_x_array = zeros(1,num_waypoints+1); % Car RF X Coordinates (Updated from 0 as code runs)
car_RR_x_array = zeros(1,num_waypoints+1); % Car RR X Coordinates (Updated from 0 as code runs)
car_LF_x_array = zeros(1,num_waypoints+1); % Car LF X Coordinates (Updated from 0 as code runs)
car_LR_x_array = zeros(1,num_waypoints+1); % Car LR X Coordinates (Updated from 0 as code runs)
car_RF_y_array = zeros(1,num_waypoints+1); % Car RF Y Coordinates (Updated from 0 as code runs)
car_RR_y_array = zeros(1,num_waypoints+1); % Car RR Y Coordinates (Updated from 0 as code runs)
car_LF_y_array = zeros(1,num_waypoints+1); % Car LF Y Coordinates (Updated from 0 as code runs)
car_LR_y_array = zeros(1,num_waypoints+1); % Car LR Y Coordinates (Updated from 0 as code runs)

for i = 0:num_waypoints
    % Front Straightaway
    if i == 0 % Initial Start at (0,0)
        x_track_array(i+1) = 0; % X Coordinate of Track Instantaneously
        y_track_array(i+1) = 0; % Y Coordinate of Track Instantaneously
        theta_track_array(i+1) = 0; % Angle of Track Curvature Instantaneously
        car_RF_x_array(i+1) = car_length; % X Coordinate RF of Car Instantaneously
        car_RR_x_array(i+1) = -car_length; % X Coordinate RR of Car Instantaneously
        car_LF_x_array(i+1) = car_length; % X Coordinate LF of Car Instantaneously
        car_LR_x_array(i+1) = -car_length; % X Coordinate LR of Car Instantaneously
        car_RF_y_array(i+1) = -car_width; % Y Coordinate RF of Car Instantaneously
        car_RR_y_array(i+1) = -car_width; % Y Coordinate RR of Car Instantaneously
        car_LF_y_array(i+1) = car_width; % Y Coordinate LF of Car Instantaneously
        car_LR_y_array(i+1) = car_width; % Y Coordinate LR of Car Instantaneously
    elseif (0 < i) && (i <= num_waypoints/4) % Points from (0,0) to (900,0)
        x_track_array(i+1) = (i)*delta_s;
        y_track_array(i+1) = 0;
        theta_track_array(i+1) = 0;
        car_RF_x_array(i+1) = (i)*delta_s + car_length;
        car_RR_x_array(i+1) = (i)*delta_s - car_length;
        car_LF_x_array(i+1) = (i)*delta_s + car_length;
        car_LR_x_array(i+1) = (i)*delta_s - car_length;
        car_RF_y_array(i+1) = -car_width;
        car_RR_y_array(i+1) = -car_width;
        car_LF_y_array(i+1) = car_width;
        car_LR_y_array(i+1) = car_width;
    % Turns 1 & 2
    elseif (num_waypoints/4 < i) && (i <= num_waypoints/2) % Points from (900,0) to (900,400)
        theta_track_array(i+1) = theta_track_array(i) + delta_theta;
        x_track_array(i+1) = x_track_array(i) + turn_radius * sin(theta_track_array(i+1));
        y_track_array(i+1) = turn_radius - turn_radius * cos(theta_track_array(i+1));
        % car_RF_x_array(i+1) = x_track_array(i) + (turn_radius + car_width) * sin(theta_track_array(i+1));
        % car_RR_x_array(i+1) = x_track_array(i) + (turn_radius + car_width) * sin(theta_track_array(i+1));
        % car_LF_x_array(i+1) = x_track_array(i) + (turn_radius - car_width) * sin(theta_track_array(i+1));
        % car_LR_x_array(i+1) = x_track_array(i) + (turn_radius - car_width) * sin(theta_track_array(i+1));
        % car_RF_y_array(i+1) = turn_radius - (turn_radius + car_width) * cos(theta_track_array(i+1));
        % car_RR_y_array(i+1) = turn_radius - (turn_radius + car_width) * cos(theta_track_array(i+1));
        % car_LF_y_array(i+1) = turn_radius - (turn_radius - car_width) * cos(theta_track_array(i+1));
        % car_LR_y_array(i+1) = turn_radius - (turn_radius - car_width) * cos(theta_track_array(i+1));
    % Back Straightaway
    elseif (num_waypoints/2 < i) && (i <= num_waypoints*(3/4)) % Points from (900,400) to (0,400)
        x_track_array(i+1) = straightaway_length - (i-60)*delta_s;
        y_track_array(i+1) = 2 * turn_radius;
        theta_track_array(i+1) = 0;
        car_RF_x_array(i+1) = straightaway_length - (i-60)*delta_s - car_length;
        car_RR_x_array(i+1) = straightaway_length - (i-60)*delta_s + car_length;
        car_LF_x_array(i+1) = straightaway_length - (i-60)*delta_s - car_length;
        car_LR_x_array(i+1) = straightaway_length - (i-60)*delta_s + car_length;
        car_RF_y_array(i+1) = 2*turn_radius + car_width;
        car_RR_y_array(i+1) = 2*turn_radius + car_width;
        car_LF_y_array(i+1) = 2*turn_radius - car_width;
        car_LR_y_array(i+1) = 2*turn_radius - car_width;
    % Turns 3 & 4
    else  % Points from (0,400) to (0,0)
        disp('')
        % theta_track_array(i) = (i - num_waypoints * (3 / 4) - 1) * delta_theta;
        % x_track_array(i) = x_track_array(num_waypoints * (3 / 4)) + turn_radius * sin(theta_track_array(i));
        % y_track_array(i) = turn_radius - turn_radius * cos(theta_track_array(i));
        % car_RF_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius + car_width) * sin(theta_track_array(i));
        % car_RR_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius + car_width) * sin(theta_track_array(i));
        % car_LF_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius - car_width) * sin(theta_track_array(i));
        % car_LR_x_array(i) = x_track_array(num_waypoints * (3 / 4)) + (turn_radius - car_width) * sin(theta_track_array(i));
        % car_RF_y_array(i) = 2 * turn_radius - (turn_radius + car_width) * cos(theta_track_array(i));
        % car_RR_y_array(i) = 2 * turn_radius - (turn_radius + car_width) * cos(theta_track_array(i));
        % car_LF_y_array(i) = 2 * turn_radius - (turn_radius - car_width) * cos(theta_track_array(i));
        % car_LR_y_array(i) = 2 * turn_radius - (turn_radius - car_width) * cos(theta_track_array(i));
    end
end

% Plotting the Track
figure;
plot(x_track_array, y_track_array, 'k', 'LineWidth', track_width);
hold on;

for i = 0:num_waypoints
    % Plotting the Car
    car_x_array = [car_RF_x_array(i+1),car_RR_x_array(i+1),car_LR_x_array(i+1),car_LF_x_array(i+1)];
    car_y_array = [car_RF_y_array(i+1),car_RR_y_array(i+1),car_LR_y_array(i+1),car_LF_y_array(i+1)];
    patch(car_x_array,car_y_array,'b','LineWidth',car_width);
    
    % Plotting the Car Path
    animated_line = animatedline(x_track_array, y_track_array,'Color','r','LineWidth',1);
end

% Set axis limits
axis equal;
axis([-100, straightaway_length + 100, -50, 2 * turn_radius + 50]);

% Add labels and title
xlabel('X (m)');
ylabel('Y (m)');
title('Race Track with Racecar');

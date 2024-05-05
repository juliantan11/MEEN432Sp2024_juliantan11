%Initialize Values
GenTrack;
P4init;

%% Run Simulation %%
out = sim("Project4Final_v2.slx", "StopTime", "3600");

%% Extract Data from Simulation Output %%
car_X = out.X.Data;
car_Y = out.Y.Data;
car_time = out.X.Time;
car_vel = out.veh_speed.Data;
SOC = out.SOC.Data; 
brake_viol = out.brake.Data; 

%% Analysis with raceStat %%
path.width = track_width; % Use track_width from your initialization code
path.l_st = straight_length; % Use straight_length from your initialization code
path.radius = curve_radius; % Use curve_radius from your initialization code

raceStats = raceStat(car_X, car_Y, car_time, path, out);

%% Display Results %%
disp(raceStats);

%%% Plotting Code %%%
figure;
x_min = -400;
x_max = 1300;

y_min = -400;
y_max = 800;

patch([x_min, x_max, x_max, x_min], [y_min, y_min, y_max, y_max], [0.5 0.8 0.5], 'EdgeColor', 'none');
hold on;

% Plot the track in grey
plot(x, y, 'Color', [0.5 0.5 0.5],'LineWidth', track_width); % Plot track
xlim([x_min, x_max]);
ylim([y_min, y_max]);
axis([-400 1300 -400 800])

% Plot the vehicle's path
plot(car_X, car_Y, 'r--', 'LineWidth', 2);

% Define the rectangular patch that represents the vehicle
rect_length = track_width * 5;  % Increase the patch length
rect_width = track_width * 2.5; % Increase the patch width
vehicle_patch = patch([0, rect_length, rect_length, 0], [-rect_width/2, -rect_width/2, rect_width/2, rect_width/2], 'k');
rotate(vehicle_patch, [0 0 1], rad2deg(theta(1)), [x(1) y(1) 0]);
vehicle_patch.EdgeColor = [0 0 0];
vehicle_patch.FaceColor = 'b';
alpha(vehicle_patch, 0.8); % Set transparency


axis equal;
title('Project 4: 1 Hour Vehicle Loop');
xlabel('X (m)');
ylabel('Y (m)');
legend('', 'Track', 'Vehicle Path');

% Animation loop for the vehicle patch following the waypoints
for i = 2:length(car_X)
    % Calculate new theta for orientation
    theta(i) = atan2(car_Y(i) - car_Y(i-1), car_X(i) - car_X(i-1));
    
    % Update vehicle position and orientation
    set(vehicle_patch, 'XData', vehicle_patch.XData + (car_X(i) - car_X(i-1)), ...
                       'YData', vehicle_patch.YData + (car_Y(i) - car_Y(i-1)));
    rotate(vehicle_patch, [0 0 1], rad2deg(theta(i) - theta(i-1)), [car_X(i) car_Y(i) 0]);
    
    drawnow; % Update the plot with new position and orientation
end

hold off; % Release the plot hold



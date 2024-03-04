% Vehicle Parameters
carData.Inertia = 1600; % kg m^2  -  Car Inertia
carData.Mass = 1000; % kg      -  Car Mass 

% Initial Conditions
carData.init.X0 = 0;        % m - Initial X Position of the Car
carData.init.Y0 = 0;        % m - Initial Y Position of the Car
carData.init.vx0 = 0.1;     % m/s - Initial Velocity in X of the Car
carData.init.vy0 = 0;       % m/s - Initial Velocity in Y of the Car
carData.init.omega0 = 0;    % rad/s - Initial Yaw Rate of the Car
carData.init.psi0 = 0;      % rad - Initial Heading of the Car

% Vehicle Tire Information 
carData.Calpha_f = 40000; % N/rad - Front Tire Coefficient (slope)
carData.Calpha_r = 40000; % N/rad - Rear Tire Coefficient (slope)
carData.Fyfmax = 40000*1/180*pi; % N - Max Front Tire Force
carData.Fyrmax = 40000*1/180*pi; % N - Max Rear Tire Force
carData.lr = 1.5; % m - Distance from CG to rear axis
carData.lf = 1.0; % m - Distance from CG to front axis
carData.radius = 0.3; % m - Radius of tires

% Track Specs
track.radius = 200; % Radius of Curves
track.width = 15; % Width of the Track
track.l_straightaways = 900; % Length of Straightaways

% Waypoint Setup
num_waypoints = 120; % Number of Waypoints for Track
delta_s = (2*track.l_straightaways)/(num_waypoints/2); % Change in Distance in Track Straightaways (30 for 120) [m]
delta_theta = (2*pi)/(num_waypoints/2); % Change in Angle in Track Turns (pi/30 for 120) [rad]

% Initialize arrays to store X, Y coordinates, and time
X = zeros(1, num_waypoints+1);
Y = zeros(1, num_waypoints+1);
time = zeros(1, num_waypoints+1);

% Simulate the vehicle going around the track
for i = 0:num_waypoints
    % Calculate X and Y coordinates
    if i <= num_waypoints/4
        X(i+1) = i * delta_s;
        Y(i+1) = 0;
    elseif i <= num_waypoints/2
        theta = (i - num_waypoints/4) * delta_theta;
        X(i+1) = track.l_straightaways + track.radius * sin(theta);
        Y(i+1) = track.radius - track.radius * cos(theta);
    elseif i <= num_waypoints*(3/4)
        X(i+1) = track.l_straightaways - (i - num_waypoints/2) * delta_s;
        Y(i+1) = 2 * track.radius;
    else
        theta = (i - num_waypoints/2) * delta_theta;
        X(i+1) = -track.radius * sin(theta);
        Y(i+1) = track.radius + track.radius * cos(theta);
    end
    
    % Calculate time (for simplicity, assuming constant speed)
    time(i+1) = i; % Adjust as needed based on your simulation
    
    % Plot the track and vehicle
    plot(X, Y, 'k', 'LineWidth', track.width);
    hold on;
    plot(X(i+1), Y(i+1), 'ro'); % Plot vehicle position
    hold off;
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    title('Vehicle Going Around the Track');
    drawnow;
end

% Run raceStat function to analyze the data
raceStats = raceStat(X, Y, time, track);

% Display results
disp('Number of loops completed: ');
disp(raceStats.loops);
disp('Completion time: ');
disp(raceStats.tloops(end)); % Assuming last element in tloops is the completion time
if isempty(raceStats.leftTrack.X)
    disp('Vehicle stayed on track.');
else
    disp('Vehicle went off track.');
end

% Plot Setup of Track and Car
xpath = zeros(1,num_waypoints+1); % X Track Coordinates (Updated from 0 as code runs)
ypath = zeros(1,num_waypoints+1); % Y Track Coordinates (Updated from 0 as code runs)
thetapath = zeros(1,num_waypoints+1); % Angle Track Coordinates (Updated from 0 as code runs)

for i = 0:num_waypoints
    % Front Straightaway
    if i == 0 % Initial Start at (0,0)
        xpath(i+1) = 0; % X Coordinate of Track Instantaneously
        ypath(i+1) = 0; % Y Coordinate of Track Instantaneously
        thetapath(i+1) = 0; % Angle of Track Curvature Instantaneously
        
    elseif (0 < i) && (i <= num_waypoints/4) % Points from (0,0) to (900,0)
        xpath(i+1) = (i)*delta_s;
        ypath(i+1) = 0;
        thetapath(i+1) = 0;
        
    % Turns 1 & 2
    elseif (num_waypoints/4 < i) && (i <= num_waypoints/2) % Points from (900,0) to (900,400)
        thetapath(i+1) = thetapath(i) + delta_theta;
        xpath(i+1) = straightaway_length + track.radius*sin(thetapath(i+1));
        ypath(i+1) = track.radius - track.radius*cos(thetapath(i+1));
        
    % Back Straightaway
    elseif (num_waypoints/2 < i) && (i <= num_waypoints*(3/4)) % Points from (900,400) to (0,400)
        xpath(i+1) = straightaway_length - (i-60)*delta_s;
        ypath(i+1) = 2 * track.radius;
        thetapath(i+1) = 0;
        
    % Turns 3 & 4
    else  % Points from (0,400) to (0,0)
        thetapath(i+1) = thetapath(i) + delta_theta;
        xpath(i+1) = -track.radius*sin(thetapath(i+1));
        ypath(i+1) = track.radius + track.radius*cos(thetapath(i+1));
        
    end
end

% Plotting the Track
figure;
plot(xpath, ypath, 'k', 'LineWidth', track.width);
hold on;

% Plotting the Car Path
carpath = animatedline('Color','r');

% Set axis limits
axis equal;
axis([-100, straightaway_length + 100, -50, 2 * track.radius + 50]);

% Add labels and title
xlabel('X (m)');
ylabel('Y (m)');
title('Race Track with Racecar');

% create a "car" of width w and length 2w
w = 0.2;
carpos = [- w/2, - w; w/2, -w; w/2, w; -w/2, w]';
a = patch('XData',carpos(:,1),'YData',carpos(:,2));
a.EdgeColor = [0 0 1];
a.FaceColor = 'blue';

% perform an animated "simulation" - no dynamics, just kinematics
for i = 0:num_waypoints
    % Plotting the Car Path
    addpoints(carpath, xpath(i+1), ypath(i+1))
    drawnow
    
    % Front Straightaway
    if i == 0 % Initial Start at (0,0)
        carpos = [- w/2, - w; w/2, -w; w/2, w; -w/2, w]'; % Plots Initial Car
        
    elseif (0 < i) && (i <= num_waypoints/4) % Points from (0,0) to (900,0)
        carpos = [xpath(i+1) - w/2, ypath(i+1) - w; xpath(i+1) + w/2, ypath(i+1) - w; xpath(i+1) + w/2, ypath(i+1) + w; xpath(i+1) - w/2, ypath(i+1) + w]; % Plot Car in Straights
        
    % Turns 1 & 2
    elseif (num_waypoints/4 < i) && (i <= num_waypoints/2) % Points from (900,0) to (900,400)
        carpos = rotate(carpos'-[xpath(i+1);ypath(i+1)], 0)' + [xpath(i+1);ypath(i+1)]'; % Plot Car in Turns

        
    % Back Straightaway
    elseif (num_waypoints/2 < i) && (i <= num_waypoints*(3/4)) % Points from (900,400) to (0,400)
        carpos = [xpath(i+1) - w/2, ypath(i+1) - w; xpath(i+1) + w/2, ypath(i+1) - w; xpath(i+1) + w/2, ypath(i+1) + w; xpath(i+1) - w/2, ypath(i+1) + w]; % Plot Car in Straights
        
    % Turns 3 & 4
    else  % Points from (0,400) to (0,0)
        carpos = rotate(carpos'-[xpath(i+1);ypath(i+1)], 0)' + [xpath(i+1);ypath(i+1)]'; % Plot Car in Turns

        
    end

    %update car image
    a.XData = carpos(:,1);
    a.YData = carpos(:,2);
    drawnow
end

hold off

% The following lines of code should be moved below the last function definition
XYposplot = sim('Project2Week2_.slx');

plot(XYposplot)

% Move this below all function definitions
function raceStats = raceStat(X,Y,t,path)
%========================================================
% 
% Usage: rs = raceStat(X,Y,t,path)
%
% Inputs:   X, Y are coordinates from your vehicle simulations. 
%           t is the set times corresponding to X and Y
%           path is a structure of with fields "width" (width of the 
%                  track), "l_st" (length of the straight away), and 
%                  "radius" (radius of the curved section)
%
% Outputs: raceStats is a structure with the following fields:
%
%   loops - this is an integer that tells how many loops around
%              the track the vehicle has gone around
%   tloops - this is an array that tells you the time(s) that the
%              start line was crossed 
%   lefftTrack - this has the X,Y and t values when the vehicle
%              went outside the track
%   
%========================================================

prev_section = 6;
loops = -1;
j = 0;
k = 0;
Xerr = [];
Yerr = [];
terr = [];
for i = 1:length(X)
    if X(i) < path.l_straightaways
        if X(i) >= 0
            if Y(i) < path.radius
                section = 1;
            else
                section = 4;
            end
        else
            if Y(i) < path.radius
                section = 6;
            else
                section = 5;
            end
        end
    else
        if Y(i) < path.radius
            section = 2;
        else
            section = 3;
        end
    end
    if ((prev_section == 6) && (section == 1))
        loops = loops  + 1;
        j = j+1;
        tloops(j) = t(i);
    end
    prev_section = section;
    if ~insideTrack(X(i),Y(i),section,path)
        k = k+1;
        Xerr(k) = X(i);
        Yerr(k) = Y(i);
        terr(k) = t(i);
    end
end
raceStats.loops = loops;
raceStats.tloops = tloops;
raceStats.leftTrack.X = Xerr;
raceStats.leftTrack.Y = Yerr;
raceStats.leftTrack.t = terr;
end

function yesorno = insideTrack(x,y,section,path)
switch section
    case 1
        if ((y < (0.0 + path.width)) && (y > (0.0 - path.width))) 
            yesorno = 1;
        else
            yesorno = 0;
        end
    case {2, 3}
        rad = sqrt((x - path.l_straightaways)^2 + (y - path.radius)^2); % Corrected field name
        if ((rad < path.radius + path.width) && ...
                (rad > path.radius - path.width))
            yesorno = 1;
        else
            yesorno = 0;
        end
    case 4
        if ((y < (2 * path.radius + path.width)) && ...
                (y > (2 * path.radius - path.width))) 
            yesorno = 1;
        else
            yesorno = 0;
        end        
    case {5, 6}
        rad = sqrt((x - 0.0)^2 + (y - path.radius)^2);
        if ((rad < path.radius + path.width) && ...
                (rad > path.radius - path.width))
            yesorno = 1;
        else
            yesorno = 0;
        end
    otherwise
        print("error");
end
end

function xyt = rotate(xy,theta)
xyt = TF(theta) * xy;
end

function y = TF(psi)
y = [cos(psi), sin(psi); -sin(psi), cos(psi)];
end

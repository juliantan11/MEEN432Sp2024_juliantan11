% Project Week 1 Demo

% Initial Conditions
w_0 = 1; % Initial Angular Velocity [rad/s]
A = 1; % Constant Applied Torque [N*m]
b = 1; % Damping Coefficient [ N*m*s/rad]
J = 1; % Rotational Inertia [kg*m^2]

dT = [0.001, 0.1, 1]; % Time Step [s]
solvers = ["ode1", "ode4"]; % Fixed Time Step Solver [Euler]

% Create figure for plotting
figure;

% Nested loops to iterate over time steps and solvers
for i = 1:length(dT)
    for j = 1:length(solvers)
        solver = solvers(j);
        dT_val = dT(i);

        % Simulate the system
        simout = sim("ProjectWeek_1.slx", "Solver", solver, "FixedStep", string(dT_val));

        % Extract data
        W = simout.w.Data;
        W_DOT = simout.w_dot.Data;
        T = simout.tout;

        % Calculate subplot index
        subplot_index = (i - 1) * length(solvers) + j;

        % Plot angular velocity
        subplot(length(dT), length(solvers)*2, subplot_index);
        plot(T, W);
        title(['Solver: ', char(solver), ', Time Step: ', num2str(dT_val)]);
        xlabel('Time [s]');
        ylabel('Angular Velocity [rad/s]');

        % Plot angular acceleration
        subplot(length(dT), length(solvers)*2, subplot_index + length(solvers)*2); 
        plot(T, W_DOT);
        title(['Solver: ', char(solver), ', Time Step: ', num2str(dT_val)]); % Added title
        xlabel('Time [s]');
        ylabel('Angular Acceleration [rad/s^2]');
    end
end

% Initial Conditions
w_0 = 1; % Initial Angular Velocity [rad/s]
A = 1; % Constant Applied Torque [N*m]
b = 1; % Damping Coefficient [ N*m*s/rad]
J = 1; % Rotational Inertia [kg*m^2]

dT = [0.001, 0.1, 1]; % Time Step [s]



% Project Week 2
solvers = ["ode1", "ode4", "ode45", "ode23tb"];

ode1_dt = {}
ode1_max_error = {}
ode1_cpu_time = {}

ode4_dt = {}
ode4_max_error = {}
ode4_cpu_time = {}

ode45_max_error = {}
ode45_cpu_time = {}

ode23tb_max_error = {}
ode23tb_cpu_time = {}

for i = 1:length(solvers)
    if solver{i,1}, Solver == "ode1" or "ode4"
        dT = [0.001, 0.1, 1]; % Fixed Time Step Values [s]
        for j = 1:length(dT)
            A = [0, 100]; % Constant Torque Values [N*m]
            for k = 1:length(A)
                b = [10, 0.1]; % Damping Coefficient [N*m*s/rad]
                for l = 1:length(b)
                    J1 = [100, 0.01]; % Rotational Inertia [kg*m^2]
                    for m = 1:length(J1)
                        w_0 = [10, 0.0]; % Initial Conditions [rad/s]
                        for n = 1:length(w_0)
                            dT_val = j;
                            A_val = k;
                            b_val = l;
                            J1_val = m;
                            w_0_val = n;
                            simout = sim("Project1Week2_.slx", "Solver", solver, "FixedStep", string(dT_val));
                            
    elseif srting(i) == "ode45" or "ode23tb"
        F = [0.1, 100]; % Frequency of Torque Values [rad/s]
        for i = 1:length(F)
            
        
        
    else
        print("Error solver selection")
    end
end




% Create figure for plotting
figure;

% Nested loops to iterate over time steps and solvers
for i = 1:length(dT)
    for j = 1:length(solvers)
        solver = solvers(j);
        dT_val = dT(i);

        % Simulate the system
        simout = sim("Project1Week2_.slx", "Solver", solver, "FixedStep", string(dT_val));

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

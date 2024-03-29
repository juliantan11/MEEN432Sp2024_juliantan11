% Part 1

solver_arr = {'ode1', 'ode4', 'ode45', 'ode23tb'};

% Initialize arrays to store simulation results
ode1_results = struct('dT', [], 'max_error', [], 'cpu_time', []);
ode4_results = struct('dT', [], 'max_error', [], 'cpu_time', []);
ode45_results = struct('max_error', [], 'cpu_time', []);
ode23tb_results = struct('max_error', [], 'cpu_time', []);

for i = 1:length(solver_arr)
    solver = solver_arr{i};

    if strcmp(solver, 'ode1') || strcmp(solver, 'ode4')
        dT_arr = [0.001, 0.1, 1]; % Fixed Time Step Values [s]
        A_arr = [0, 100]; % Constant Torque Values [N*m]
        b_arr = [10, 0.1]; % Damping Coefficient [N*m*s/rad]
        J1_arr = [100, 0.01]; % Rotational Inertia [kg*m^2]
        w_0_arr = [10, 0.0]; % Initial Conditions [rad/s]
        isSin = 0; % Boolean Operator for Sine Function Initialization [False]
        F = 0; % Frequency [rad/s]

        for dT = dT_arr
            for A = A_arr
                for b = b_arr
                    for J1 = J1_arr
                        for w_0 = w_0_arr
                            % Start CPU time
                            cpu_time_start = cputime;

                            % Simulate the system
                            simout = sim('Project1Part1_.slx', ...
                                'Solver', solver, 'FixedStep', num2str(dT));

                            % End CPU time
                            cpu_time_end = cputime;
                            cpu_time = cpu_time_end - cpu_time_start;

                            % Extract data
                            W = simout.w.Data;
                            T = simout.tout;

                            % Error Calculation
                            W_len = length(W);
                            dt = linspace(0, 25, W_len);
                            theory_w = theory_omega(dt, A, b, w_0, J1, isSin);
                            theory_w = transpose(theory_w);
                            max_error = max(abs(W - theory_w));

                            % Store results
                            if strcmp(solver, 'ode1')
                                ode1_results.dT(end+1) = dT;
                                ode1_results.cpu_time(end+1) = cpu_time;
                                ode1_results.max_error(end+1) = max_error;
                            else
                                ode4_results.dT(end+1) = dT;
                                ode4_results.cpu_time(end+1) = cpu_time;
                                ode4_results.max_error(end+1) = max_error;
                            end
                        end
                    end
                end
            end
        end
    elseif strcmp(solver, 'ode45') || strcmp(solver, 'ode23tb')
        A_arr = [0, 100]; % Constant Torque Values [N*m]
        b_arr = [10, 0.1]; % Damping Coefficient [N*m*s/rad]
        J1_arr = [100, 0.01]; % Rotational Inertia [kg*m^2]
        w_0_arr = [10, 0.0]; % Initial Conditions [rad/s]
        isSin = 1; % Boolean Operator for Sine Function Initialization [True]
        F_arr = [0.1, 100]; % Frequency [rad/s]

        for A = A_arr
            for b = b_arr
                for J1 = J1_arr
                    for w_0 = w_0_arr
                        for F = F_arr
                            % Start CPU time
                            cpu_time_start = cputime;

                            % Simulate the system
                            simout = sim('Project1Part1_.slx', ...
                                'Solver', solver);

                            % End CPU time
                            cpu_time_end = cputime;
                            cpu_time = cpu_time_end - cpu_time_start;

                            % Extract data
                            W = simout.w.Data;
                            T = simout.tout;

                            % Error Calculation
                            W_len = length(W);
                            dt = linspace(0, 25, W_len);
                            theory_w = theory_omega(dt, A, b, w_0, J1, isSin);
                            theory_w = transpose(theory_w);
                            max_error = max(abs(W - theory_w));

                            % Store results
                            if strcmp(solver, 'ode45')
                                ode45_results.cpu_time(end+1) = cpu_time;
                                ode45_results.max_error(end+1) = max_error;
                            else
                                ode23tb_results.cpu_time(end+1) = cpu_time;
                                ode23tb_results.max_error(end+1) = max_error;
                            end
                        end
                    end
                end
            end
        end
    else
        disp(['Solver ', solver, ' not implemented']);
    end
end

% Display results
disp('ode1 results:');
disp(ode1_results);
disp('ode4 results:');
disp(ode4_results);
disp('ode45 results:');
disp(ode45_results);
disp('ode23tb results:');
disp(ode23tb_results);

% Plot 1: Max Error vs Time Step (Fixed Time Step)
figure;
plot(ode1_results.dT, ode1_results.max_error, '-o', 'DisplayName', 'ODE1');
hold on;
plot(ode4_results.dT, ode4_results.max_error, '-s', 'DisplayName', 'ODE4');
xlabel('Time Step [s]');
ylabel('Max Error');
title('Max Error vs Time Step (Fixed Time Step Solvers)');
legend('show');
grid on;

% Plot 2: CPU Time vs Time Step (Fixed Time Step)
figure;
plot(ode1_results.dT, ode1_results.cpu_time, '-o', 'DisplayName', 'ODE1');
hold on;
plot(ode4_results.dT, ode4_results.cpu_time, '-s', 'DisplayName', 'ODE4');
xlabel('Time Step [s]');
ylabel('CPU Time Taken [s]');
title('CPU Time vs Time Step (Fixed Time Step Solvers)');
legend('show');
grid on;

% Plot 3: Max Simulation Error vs CPU Time Taken for Fixed Time Step solvers
figure;
scatter(ode1_results.cpu_time, ode1_results.max_error, 'o', 'DisplayName', 'ODE1');
hold on;
scatter(ode4_results.cpu_time, ode4_results.max_error, 's', 'DisplayName', 'ODE4');
hold on;
scatter(ode45_results.cpu_time, ode45_results.max_error, 'd', 'DisplayName', 'ODE45');
hold on;
scatter(ode23tb_results.cpu_time, ode23tb_results.max_error, '^', 'DisplayName', 'ODE23TB');
xlabel('CPU Time Taken [s]');
ylabel('Max Error');
title('Max Simulation Error vs CPU Time Taken (Fixed Time Step Solvers)');
legend('show');
grid on;


% Part 2

solver2_arr = {'ode1', 'ode4', 'ode45'};

% Initialize arrays to store simulation results
ode1_results = struct('omega_shaft', [], 'cpu_time', []);
ode4_results = struct('omega_shaft', [], 'cpu_time', []);
ode45_results = struct('omega_shaft', [], 'cpu_time', []);

for i = 1:length(solver_arr)
    solver = solver_arr{i};
    
    if strcmp(solver, 'ode1') || strcmp(solver, 'ode4')
        dT_arr = [0.1, 1]; % Fixed Time Step Values [s]
        A_arr = [1, 100]; % Constant Torque Values [N*m]
        b1_arr = 1; % Damping Coefficient [N*m*s/rad]
        J1_arr = 100; % Rotational Inertia [kg*m^2]
        b2_arr = 1; % Damping Coefficient [N*m*s/rad]
        J2_arr = 1; % Rotational Inertia [kg*m^2]
        k_arr = [10, 100, 1000]; % Rotational Inertia [kg*m^2]
        w_0_arr = [10, 0.0]; % Initial Conditions [rad/s]

        for dT = dT_arr
            for A = A_arr
                for b1 = b1_arr
                    for J1 = J1_arr
                        for b2 = b2_arr
                            for J2 = J2_arr
                                for k = k_arr
                                    for w_0 = w_0_arr
                                        % Start CPU time
                                        cpu_time_start = cputime;
    
                                        % Simulate the system
                                        simout = sim('Project1Part2_.slx', ...
                                            'Solver', solver, 'FixedStep', num2str(dT));
    
                                        % End CPU time
                                        cpu_time_end = cputime;
                                        cpu_time = cpu_time_end - cpu_time_start;
    
                                        % Extract data
                                        out_total = simout.out_total.Data;
                                        T = simout.tout;
    
                                        % Store results
                                        if strcmp(solver, 'ode1')
                                            ode1_results.cpu_time(end+1) = cpu_time;
                                        else
                                            ode4_results.cpu_time(end+1) = cpu_time;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif strcmp(solver, 'ode45')
        A_arr = [1, 100]; % Constant Torque Values [N*m]
        b1_arr = 1; % Damping Coefficient [N*m*s/rad]
        J1_arr = 100; % Rotational Inertia [kg*m^2]
        b2_arr = 1; % Damping Coefficient [N*m*s/rad]
        J2_arr = 1; % Rotational Inertia [kg*m^2]
        k_arr = [10, 100, 1000]; % Rotational Inertia [kg*m^2]
        w_0_arr = [10, 0.0]; % Initial Conditions [rad/s]
        
        for A = A_arr
            for b = b_arr
                for J1 = J1_arr
                    for k = k_arr
                        for w_0 = w_0_arr
                            % Start CPU time
                            cpu_time_start = cputime;
    
                            % Simulate the system
                            simout = sim('Project1Part2_.slx', ...
                                'Solver', solver);
    
                            % End CPU time
                            cpu_time_end = cputime;
                            cpu_time = cpu_time_end - cpu_time_start;
    
                            % Extract data
                            out_total = simout.out_total.Data;
                            T = simout.tout;
    
                            % Store results
                            ode45_results.cpu_time(end+1) = cpu_time;
                        end
                    end
                end
            end
        end
    else
        disp(['Solver ', solver, ' not implemented']);
    end
end

figure;
plot(out_total,cpu_time, '-o', 'DisplayName', 'k=10')
hold on;
plot(out_total,cpu_time, 's', 'DiplayName', 'k=100' )
plot(out_total, cpu_time, 'd', 'DisplayName', 'k=1000')
xlabel('Time')
ylabel('Shaft Speed')
title('Shaft Speed vs. Time')
legend('show')
grid on;

figure;
plot(out_total,cpu_time, '-o', 'DisplayName', 'torque=1')
hold on;
plot(out_total,cpu_time, 's', 'DiplayName', 'torque=100' )
xlabel('Time')
ylabel('Shaft Speed')
title('Shaft Speed vs. Time')
legend('show')
grid on;

% Theory Omega
function theory_w = theory_omega(dt, tau, b, w_0, J1, isSin)
    if isSin == 0
        theory_w = (tau/b)*(1 - exp(-b*dt/J1)) + w_0*exp(-b*dt/J1);
    else
        % % Simulate the system
        solver = 'ode4';
        dT = 0.001;
        simout = sim('Project1Part1_.slx', ...
            'Solver', solver, 'FixedStep', num2str(dT));

        % Extract data
        W_ode4 = simout.w.Data;

        W_ode4 = (tau/b)*(1 - exp(-b*dt/J1)) + w_0*exp(-b*dt/J1);
        theory_w = W_ode4;
    end
end

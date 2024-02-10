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
        
        for dT_val = dT_arr
            for A = A_arr
                for b = b_arr
                    for J1 = J1_arr
                        for w_0 = w_0_arr
                            % Start CPU time
                            cpu_time_start = cputime;

                            % Simulate the system
                            simout = sim('Project1Week2_.slx', ...
                                'Solver', solver, 'FixedStep', num2str(dT_val));

                            % End CPU time
                            cpu_time_end = cputime;
                            cpu_time = cpu_time_end - cpu_time_start;

                            % Extract data
                            W = simout.w.Data;
                            T = simout.tout;

                            % Error Calculation
                            W_len = length(W)
                            dt = linspace(0, 25, W_len)
                            theory_w = theory_w(dt, A, b, w_0, J1)
                            error = cat(1, error, max(W - theory_w));
                            
                            % Plot angular velocity
                            %figure;
                            %plot(T, W);
                            %title(['Solver: ', solver, ', Time Step: ', num2str(dT_val)]);
                            %xlabel('Time [s]');
                            %ylabel('Angular Velocity [rad/s]');
                            
                            % Store results
                            if strcmp(solver, 'ode1')
                                ode1_results.dT(end+1) = dT_val;
                                ode1_results.cpu_time(end+1) = cpu_time;
                            else
                                ode4_results.dT(end+1) = dT_val;
                                ode4_results.cpu_time(end+1) = cpu_time;
                            end
                        end
                    end
                end
            end
        end
    elseif strcmp(solver, 'ode45')
        % Code for ode45 solver
    elseif strcmp(solver, 'ode23tb')
        % Code for ode23tb solver
    else
        disp(['Solver ', solver, ' not implemented']);
    end
end

% Theory Omega
function theory_w(dt, tau, b, w_0, J1)
    theory_w = (tau/b)*(1 - exp(-b*T/J1)) + w_0*exp(-b*T/J1);
end

% Display results
disp('ode1 results:');
disp(ode1_results);
disp('ode4 results:');
disp(ode4_results);

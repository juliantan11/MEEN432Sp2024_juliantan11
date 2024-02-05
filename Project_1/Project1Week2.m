% Project Week 2
solver_arr = {"ode1", "ode4", "ode45", "ode23tb"};

% Fixed Time Step Arrays
ode1_dT = {};
ode1_max_error = {};
ode1_cpu_time = {};

ode4_dT = {};
ode4_max_error = {};
ode4_cpu_time = {};

% Variable Time Step Arrays
ode45_max_error = {};
ode45_cpu_time = {};

ode23tb_max_error = {};
ode23tb_cpu_time = {};

for i = 1:length(solver_arr)
    if solver_arr{i,1} == "ode1"
        dT_arr = [0.001, 0.1, 1]; % Fixed Time Step Values [s]
        A_arr = [0, 100]; % Constant Torque Values [N*m]
        b_arr = [10, 0.1]; % Damping Coefficient [N*m*s/rad]
        J1_arr = [100, 0.01]; % Rotational Inertia [kg*m^2]
        w_0_arr = [10, 0.0]; % Initial Conditions [rad/s]
        for j = 1:length(dT_arr)
            for k = 1:length(A_arr)
                for l = 1:length(b_arr)
                    for m = 1:length(J1_arr)
                        for n = 1:length(w_0_arr)
                            % Initial Conditions
                            solver = solver_arr{i,1};
                            dT = dT_arr(j);
                            A = A_arr(k);
                            b = b_arr(l);
                            J1 = J1_arr(m);
                            F = 0;
                            w_0 = w_0_arr(n);
                            
                            % Simulate the system
                            simout = sim("Project1Week2_.slx", "Solver", solver, "FixedStep", string(dT));

                            % Extract data
                            W = simout.w.Data;
                            W_DOT = simout.w_dot.Data;
                            T = simout.tout;
                            dT_simout = simout(i,1).SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize;
                            
                            % Error Calculation
                            % theory_w = 
                            % w_error = theory_w - W
                            
                            % CPU Time
                            % 

                            ode1_dT = cat(1, ode1_dT, dT_simout);
                            % ode1_max_error.append(w_error)
                            % ode1_cpu_time.append()
                            
                            % Calculate subplot index
                            %subplot_index = (i - 1) * length(solvers) + j;
                            
                            % Plot angular velocity
                            %subplot(length(dT), length(solvers)*2, subplot_index);
                            %plot(T, W);
                            %title(['Solver: ', char(solver), ', Time Step: ', num2str(dT_val)]);
                            %xlabel('Time [s]');
                            %ylabel('Angular Velocity [rad/s]');
                            
                            % Plot angular acceleration
                            %subplot(length(dT), length(solvers)*2, subplot_index + length(solvers)*2); 
                            %plot(T, W_DOT);
                            %title(['Solver: ', char(solver), ', Time Step: ', num2str(dT_val)]); % Added title
                            %xlabel('Time [s]');
                            %ylabel('Angular Acceleration [rad/s^2]');
                        end
                    end
                end
            end
        end
    
    % elseif string(i) == "ode45" or "ode23tb"
    %     F = [0.1, 100]; % Frequency of Torque Values [rad/s]
    %     for i = 1:length(F)

    %     end
    else
        
    end
end

print(ode1_dT)

% Project Week 2
solver = ["ode1", "ode4", "ode45", "ode23tb"];

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

for i = 1:length(solver)
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
                            % Initial Conditions
                            solver_val = solver(i);
                            dT_val = dT(j);
                            A_val = A(k);
                            b_val = b(l);
                            J1_val = J1(m);
                            w_0_val = w_0(n);
                            
                            % Simulate the system
                            simout = sim("Project1Week2_.slx", "Solver", solver_val, "FixedStep", string(dT_val));

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
                end
            end
        end
    
    elseif srting(i) == "ode45" or "ode23tb"
        F = [0.1, 100]; % Frequency of Torque Values [rad/s]
        for i = 1:length(F)
            
        
        
    else
        print("Error solver selection")
    end
end

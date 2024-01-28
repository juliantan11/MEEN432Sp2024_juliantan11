%Project Week 1 Demo

%Initial Conditions
w_0 = 1; % Initial Angular Velocity [rad/s]
A = 4; % Constant Applied Torque [N*m]
b = 1; % Damping Coefficient [ N*m*s/rad]
J = 1; % Rotational Inertia [kg*m^2]
dT = [0.001,0.1,1]; % Time Step [s]
solver = ["ode1", "ode4"]; % Fixed Time Step Solver [Euler]
simout = sim("ProjectWeek_1.slx", "Solver", solver, "FixedStep", string(dT)); 

W = simout.w.Data;
W_DOT = simout.w_dot.Data;
T = simout.tout;

plot(W,T); 
plot(W_DOT,T);
% Project 3 Final Retry

% Urban Drive Cycle
DriveData = DriveData_Urban;
Time = Time_Urban;
Time_Max = 1369; %[s]

simout = sim("Project3_Final_Retry_.slx", "StopTime", num2str(Time_Max));
sim_vel_mps = simout.vel.Data; %[mps]
sim_vel_mph = sim_vel_mps/mph2mps; %[mph]
sim_time = simout.tout; %[s] 

figure
plot(sim_time, sim_vel_mph, 'b')
hold on
plot(Time, DriveData, '--r')
xlabel("Time (s)")
ylabel("Velocity (mph)")
legend("Sim Velocity", "Drive Cycle Velocity")
title("Sim Vehicle Velocity v. Time for Urban Drive Cycle")

% Highway Drive Cycle
DriveData = DriveData_Highway;
Time = Time_Highway;
Time_Max = 765; %[s]

simout = sim("Project3_Final_Rerun.slx", "StopTime", num2str(Time_Max));
sim_vel_mps = simout.vel.Data; %[mps]
sim_vel_mph = sim_vel_mps/mph2mps; %[mph]
sim_time = simout.tout; %[s] 

figure
plot(sim_time, sim_vel_mph, 'b')
hold on
plot(Time, DriveData, '--r')
xlabel("Time (s)")
ylabel("Velocity (mph)")
legend("Sim Velocity", "Drive Cycle Velocity")
title("Sim Vehicle Velocity v. Time for Highway Drive Cycle")

disp(MotorPower)
disp(MotorEnergy)

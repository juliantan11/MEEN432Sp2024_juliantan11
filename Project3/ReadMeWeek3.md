Utilizing the model we created for week 1 and 2, the team created a module for the electric motor drive module which contains the electric motor and battery as its subsystems. Within the electric motor subsystem, this system calls for throttle command and motor speed as its inputs and outputs motor torque and motor inertia. The battery subsystem takes current as its input into the system and outputs battery voltage. Similar to week 1 and week 2, to run this code you must first run the initialization code and run whatever drive cycle is desired. You can verify that the drive cycle selected ran correctly by clicking on the output graphs. 
For our regenerative braking system, we decided to create a subsystem within our powertrain scope to use regenerative braking. Essentially when torque is negative, our power inverter (regenerative braking subsystem) will send the current
back into the battery pack. A gain block was added to adjust for the energy lost during the regen braking.

We did not adjust the brake pedal as it is dependent on velocity of the car. Our powertrain system accounts for the negative torque from regen braking, thus allowing it to slow down the car when it is not accelerating. 

## Final Submission Feedback (39/85)
Regarding the README.md file:
1) You need to provide proper instructions for how you want me to run the project files. Saying to just run the initialization files and then run the desired drive cycle is not enough since there are multiple files that have initialization information. Also due to the resubmission there are two different simulink models with the name "Final" in them so for future submissions please be more specific.

Regarding the scripts:
1) There is no script that plots the simulated vehicle velocity of both drive cycles against each of the drive cycles showing whether or not the vehicle can stay within a +- 3 mph error band

In this feedback I will be going through each of the main components of the Simulink model and provide any corrections that should be made.

First off, the model StopTime needs to be changed for both drive cycles since they each have a different stop time of 765 and 1369 (highway and urban) so you need to either tell me in the instructions to change the stop time depending on what drive cycle is being used, or you can change the stop time using a script
1) Drive Schedule: No comments
2) Driver Model: I noticed that the regen braking logic is being done in the Powertrain system but the logic we want to use is fairly straightforward and we want to have it within the driver model.
- Essentially for regen braking, you want to take a small percentage (~10%) of your total brake%cmd and add it to the throttle%cmd, then you pass the rest of the brake% as friction brakes (FBPP(%))
3) Battery: Other than the regen braking poriton, there is also some errors with the logic of calculating the SOC.
- Recall that the equation for SOC is as follows: SOC = 1 - integral(Icell)/C, where C should be the capacity per cell in units of Amp-s. The C value given in the initialization file is for the total battery capacity in Amp-hr, so make sure to convert this total battery capacity by multipying by (3600/numofcellsinParallel)
4) Electric Motor: There are many issues with the motor
- For starters, you cannot use a constant block and call out.efficiency, or out.maxtorquewV since you are attempting to call simulation output data that does not exist when the simulation is currently running
- A lot of calculations are separated into different matlab function blocks which is not necessary as they can all be combined into a single matlab function (Calculation of Pm, sgn(Pm), Pe, etc)
- Also you have the input of motor speed, but you don't actually have it connected to anything. Additionally, when using the motor speed for the lookup tables you want to make sure that the motor speed is converted from rad/s to rpm since the motor data given in the initialization file are for motor speeds in units of rpm.
5) Transmission: There are some things wrong with the logic here
- You have one part of the transmission correct which is multiplying the motor torque by some gear ratio (either FDG or others that were given) and that is the torque going to the wheels
- Whats missing from this subsystem is the input of the wheel speed (Ww) and multipying that by the inverse of some gear ratio, which will give you the motor speed that should going into the Electric motor drive system
6) Wheel Model: It seems that this system was not updated since Inertia of the wheel (Iw) and radius of the wheel (r) are not being called correctly from the carData structure.
7) Vehicle Model: Same as the wheel model, the constants used are not being called correctly
8) Brake Model: There are a few things missing from the brake logic that I want to mention real quick. The brake has two states Locked and Unlocked and we want to try to replicate that here in the brake logic
- The brake is considered to be in the Locked state when the brake%cmd == 0, so Tb = -Tw (or Tw, depends on the sign of brake%cmd)
- The brake is considered to be in the Unlocked state when brake%cmd does not == 0, so Tb = - brake%cmd * Nb,max * (Ww/(abs(Ww) + 0.001)), where Nbmax is the calibratable gain value 10000, and Ww is the angular velocity of the wheel (might be positive, again depends on sign of brake%cmd)

There are many mistakes in this model that need to fixed in order to have a working model that can used for Project 4. I understand it might be a shock to see a low grade for this project, but just know that I will allow resubmission for this Project (due 4/19 by 11:59 PM) so that the team can get a better grade. I also urge this team to reach out to me to try to set up a meeting time where we can go over what fixes need to be made so that each student can be successful for Project 4

*Note about Resubmission:* I am allowing a resubmission of this project, however depending on when this feedback gets sent out, I may not have enough time to get it regraded before Q-drop date. So please keep that in mind.

This week, the group began by debugging our week 2 code to work properly. The main issue found was turns 3 and 4 causing issues and we debugged until the code ran correctly. When moving to week 3, we created code that relays the driver's telemetry data as the car goes around the track. This data relays information like the time to complete each lap, how many laps were completed, and if the car ever goes off the track. After debugging week 2, we were able to get the track to work properly, but for week 3, the completed loops were not being output properly. Hunter worked on the coding, Julien worked on the Simulink model, and Cameron worked with both and helped where possible and completed the read me.

## Final Week Feedback (84/85)
First thing I noticed is that the raceStat function was not properly used. All you had to give it as inputs is the entire X, Y, time array from the simulation and the track structure. I took the raceStat function out of the for loop and ran it again saw that the vehicle completed 2 loops and never left the track. Another issue is with the patch simulation itself, the vehicle rectangular/square patch can barely be seen, and when it rotates on the curved sections there is a weird floating line animation. Other than that the model looks good! For the individual members to succeed in Project 4, I would advise each member to work on the driver model from Project 2 and try utilizing some path following controller (pure pursuit). A more sophisticated model will allow for the vehicle to move faster depending on where it is on the track which will be important for Project 4, as each member will need to combine this model with the Project 3 model to create a vehicle that can move around the track as fast as possible while simulating the lateral and longitudinal dynamics of the vehicle. 
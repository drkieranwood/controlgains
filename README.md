This repository is acting as a quick file transfer between a windows PC and a Unix PC.
The files within the repository are text files relating to state-space controllers for 
the AR.Drone.2.0 MAV.

There are four controller text files. Each one contains the sizes of the state-space 
system and the relevant cell values. They represent control of the X,Y,Z,W axes (where
W is the Yaw).

The idea is that Matlab on the Windows PC is used to create the control and then writes
the values to file. Then the repository is uploaded to github and then downloaded onto
the Unix PC. The state-space controller implementation in ROS can read these files and
run the AR.Drone.2.0 feedback control.

Also included in this repository is the Matlab script to output a discrete-time 
state-space system in the correct format.

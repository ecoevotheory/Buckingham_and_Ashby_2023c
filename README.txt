This ReadMe file contains information on files associated with the paper: "Separation of evolutionary timescales in coevolving species" by Buckingham & Ashby.

All files are provided "as is", without any express or implied warranty. The authors accept no liability for damages of any kind. 

Author: Lydia Buckingham, University of Bath
Contact: ljb74@bath.ac.uk
Date: 13/07/23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMMENTS

MEX files written in C# must be compiled before use, using " mex codename.c ".

Figure 1 shows phase planes which were generated manually to demonstrate a range of possible outcomes. These phase planes are not based on any specific model and so no code is provided for this figure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DESCRIPTION OF FILES

Figure 2:
"example1_simulation.m"		 	- Source code for plotting coevolutionary trajectories for the model described in Example 1	- written in matlab (R2019b). 
"example1_simulation_function.m" 	- Function which runs an evolutionary simulation for a fixed number of timesteps		- written in matlab (R2019b).
"example1_eco_dynamics_function.c"	- Function which runs ecological dynamics for a single evolutionary timestep			- written as a MEX file in C#.


Figure 3:
"example2_simulation.m"		 	- Source code for plotting coevolutionary trajectories for the model described in Example 2	- written in matlab (R2019b). 
"example2_simulation_function.m" 	- Function which runs an evolutionary simulation for a fixed number of timesteps		- written in matlab (R2019b).
"example2_eco_dynamics_function.c"	- Function which runs ecological dynamics for a single evolutionary timestep			- written as a MEX file in C#.


See code for full description and instructions for use. 
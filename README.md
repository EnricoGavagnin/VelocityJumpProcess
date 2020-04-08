
#———————————————————————————
The MATLAB codes on this repository accompany the paper:
Gavagnin, E. and Yates, C.A., 2018. Modeling persistence of motion in a crowded environment: The diffusive limit of excluding velocity-jump processes. Physical Review E, 97(3), p.032416.

#———————————————————————————
INSTRUCTIONS:

All the simulations of the agent-based models (ABMs) and the numerical solver of the partial differential equations (PDEs) are controlled by the same script: MAIN.m 


The script MAIN.m automatically recall the following functions:

ABM.m    		          : simulator of the agent-based model; 
PDE_solver.m	      	: numerical solver of the Partial differential equation;
Initial_Conditions.m 	: function which generate the initial configuration for both the ABM and the PDE;
PLOT.m			          : auxiliary function which displays the comparison between the averaged densities of ABM and PDE.

#———————————————————————————

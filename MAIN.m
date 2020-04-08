%Code to accompany the paper:
%"Modelling persistence of motion in a crowded environment: the diffusive
%limit of excluding velocity-jump processes"
%by Enrico Gavagnin and Christian A. Yates

%Created 06/10/2017
%By Enrico Gavagnin
%email: e.gavagnin@bath.ac.uk
%%

%This script is desgined to set all the parameters of the model of persistence 
%and then to launch the corresponding functions for the simulation of the 
%agent-based model (ABM) and the numerial solver of the Partial differential 
%equation (PDE). 
%%

%Save the data.
SAVE_DATA=1;

%Seleect type of plot
plot_type=2;

% Legend of plot types
% 0: no plot
% 1: Column-averaged density profiles of ABM and PDE for total population
% 2: Column-averaged density profiles of ABM and PDE for right-polarised agents

%Save the plot
SAVE_PLOT=0;

%% Parameters initialisation 

%Define the xsize and ysize
x_size=100;
y_size=100;

%Define the inital density
d=0.5;

%Define the final time 
T_final=20;

%Define the motility rate
P_m=1;

%Define the persistence parameter
phi=0.8;

%Define the reorienting rate
P_r=0.2;

%Define the velocity
v=2; 

%Type of interaction
type=3;

%Define the total number of realisations of the ABM
M=100;


%% Data/Video name

%Give a base video name
NAME='persistence_model';

%Update the video name
NAME=[NAME,'_Pm_',num2str(P_m)];

%Update the video name
NAME=[NAME,'_rep_',num2str(M)];

%Update the video name
NAME=[NAME,'_phi_',num2str(phi)];

%Update the video name
NAME=[NAME,'_v_',num2str(v)];

%Update the video name
NAME=[NAME,'_Pt_',num2str(P_r)];

%Update the video name
NAME=[NAME,'_T_',num2str(T_final)];

%Update the video name
NAME=[NAME,datestr(now,'_dd_mm_yyyy_HH_MM_SS')];

%% Initial condition
%Set the type of initial condition (IC):
IC=1;

%Legend of types of initial conditions
% 1: Uniforly polarised distribution
% 2: Only right-polarised distribution

%Use the Initial_Conditions function to set the IC of both the ABM (L_IC)
%and PDE (u_IC)
[L_IC,u_IC]=Initial_Conditions(IC,x_size,y_size,d);


%% Simulations Aagent-Based Model

%Recall the function ABM_simulator to simulate the ABM
[Rx,Lx,Ux,Dx]=ABM(type,x_size,y_size,IC,P_m,phi,P_r,v,M,T_final,d);

%% Numerical solutions PDE

%Recall the function PDE_solver to solve the system of PDEs for the
%column-averaged densities
[Rn, Ln, Un, Dn]=PDE_solver(type,x_size,u_IC,T_final,P_m,phi,P_r,v);


%% Save Data

%If the data saving is turned on, save the data
if SAVE_DATA
    %Save all the variable as data
    save([NAME,'.mat']);
end

%% PLOT
%if the data saving is turned on, save the data
if plot_type~=0
    PLOT(plot_type,SAVE_PLOT,NAME,Rx,Lx,Ux,Dx,Rn,Ln,Un,Dn,u_IC,M,phi,P_r,P_m,T_final,v)
end

%Code to accompany the paper:
%"Modelling persistence of motion in a crowded environment: the diffusive
%limit of excluding velocity-jump processes"
%by Enrico Gavagnin and Christian A. Yates

%Created 06/10/2017
%By Enrico Gavagnin
%email: e.gavagnin@bath.ac.uk
%%
% This function is desinged to set the initial condition both for the IBM
% and for the numerical solution of the PDE.

function [ L_IC,u_IC ] = Initial_Conditions(IC,x_size,y_size,seed_dens)

%% INPUT
% IC: parameter which indicates the type of initial conditions to be
%     applied
% x_size/y_size: range of the domain
% seed_dens:inital density
%
%
% IC=1: Uniform central column initial condition with uniformly distributed
% polarisations;
% IC=2: Uniform central column initial condition with only right-polarised agents;
%
%% OUTPUT: 
% L_IC: Lattice domain intialised with entrance:
%                        2: Right Mover
%                       -2: Left Mover
%                      0.5: Up Mover
%                     -0.5: Down MOver
% u_IC: Vector for the numerical solution:
%
%           u_IC=[Right Movers;Left Movers;Other Movers]
%%


%Define size of the initially populated region
seedx=floor(x_size*0.2);
seedy=y_size;

%Intial condition for the ABM
L_IC=zeros(y_size,x_size);
L_IC(y_size/2-seedy/2+1:y_size/2+seedy/2,x_size/2-seedx/2+1:x_size/2+seedx/2)=...
    floor(rand(seedy,seedx)+seed_dens*ones(seedy,seedx));

%Define the spatial discretisation step for the numerical solution
delta_x=1;

%Initial condition for the PDE
u_IC_one_pol=zeros(1,x_size/delta_x);
u_IC_one_pol((x_size/2-seedx/2)/delta_x+1:(x_size/2+seedx/2)/delta_x)=...
    seed_dens/4*ones(1,seedx/delta_x);

if IC==1
    %% Uniforly polarised distribution
    
    %Assaign the polarisations unifomrly at random for the ABM
    L_IC=L_IC.*ceil(rand(y_size,x_size)+0.5*ones(y_size,x_size));
    L_IC(L_IC==2) = -1;
    L_IC=L_IC.*ceil(rand(y_size,x_size)+0.5*ones(y_size,x_size));
    L_IC(L_IC==1) = 0.5;
    L_IC(L_IC==-1) = -0.5;
    
    %Distribute the polarisations unifomrly in the PDE
    u_IC=[u_IC_one_pol, u_IC_one_pol, u_IC_one_pol, u_IC_one_pol];
    
elseif IC==2
    %% Only right-polarised distribution
    
    %Assaign the rright-polarisation for the ABM
    L_IC=2*L_IC;
    
    %Assaign the rright-polarisation for the PDE
    u_IC=[u_IC_one_pol, zeros(size(u_IC_one_pol)), zeros(size(u_IC_one_pol)), zeros(size(u_IC_one_pol))];
    
end
end


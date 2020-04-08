%Code to accompany the paper:
%"Modelling persistence of motion in a crowded environment: the diffusive
%limit of excluding velocity-jump processes"
%by Enrico Gavagnin and Christian A. Yates

%Created 06/10/2017
%By Enrico Gavagnin
%email: e.gavagnin@bath.ac.uk
%%
%This function is designed to visualise the results of the ABM and the PDE
%and to save the plot if required.

function [] = PLOT(plot_type,SAVE_PLOT,NAME,Rx,Lx,Ux,Dx,Rn,Ln,Un,Dn,u_IC,M,phi,P_r,P_m,T_final,v)


%% INPUT
%plot_type: specify the type of plot requested (see following legend)
%SAVE_PLOT: 0 - the plot is not saved. 1 - the plot is saved with the name "NAME"
%Rx,Lx,Ux,Dx: matrices of the repeats-averaged densities of the four subpopulations of the ABM
%Rn,Ln,Dn,Un: vectors of the column-averaged densities of the four subpopulations of the PDE
%u_IC: initial condition of the numerical solution of the PDE
%M: number of ABM realisations
%phi: persistence parameter
%P_r: reorienting rate
%P_m: motility rate
%%

% Define the spatial discretisation step of the numerical solution
delta_x=1*10^(0);

%Read the size of the domain from the dimensions of the matrix Rx
[y_size, x_size]=size(Rx)

%Define the x-axises for the ABM and PDE
xn=delta_x/2:delta_x:x_size;
xx=0.5:x_size-0.1+0.5;

%Define the figure 
H=subplot(1,1,1);

if plot_type==1
    %% Column-averaged density profiles of ABM and PDE for total population
    plot(xn,(Rn+Ln+Un+Dn),'Color',[0.9 0 0],'Linewidth',3.7)
    hold on
    plot(xx,sum(Rx+Lx+Ux+Dx)/y_size,'Color',[0.1 0.1 0.1],'Linewidth',1.1)
    hold off
    title(sprintf('Comparison of the column-averaged densities of ABM (over %d repeats) and PDE \nwith P_m=%1.1f, phi=%1.1f, P_t=%1.2f, v=%d, delta=%1.1f at time %d, y-size: %d ',M,P_m,phi,P_r,v,delta_x,T_final,y_size))
    axis([0 x_size 0 4*max(sum(Rx))/y_size*1.1])
    legend('PDE','ABM')
    xlabel('Column')
    ylabel('Density')
    grid on
    
elseif plot_type==2
    %% Column-averaged density profiles of ABM and PDE for right-polarised agents
    plot(xn,Rn,'Color',[0.9 0 0],'Linewidth',3.7)
    hold on
    plot(xx,sum(Rx)/y_size,'Color',[0.1 0.1 0.1],'Linewidth',1.1)
    hold off
    title(sprintf('Comparison of the column-averaged densities of ABM (over %d repeats) and PDE \nwith P_m=%1.1f, phi=%1.1f, P_t=%1.2f, v=%d, delta=%1.1f at time %d, y-size: %d ',M,P_m,phi,P_r,v,delta_x,T_final,y_size))
    axis([0 x_size 0 max(sum(Rx))/y_size*1.1])
    legend('PDE','ABM')
    xlabel('Column')
    ylabel('Density')
    grid on
    
end

%% Save plot
if SAVE_PLOT
    SAVE_NAME=['plots/',NAME,datestr(now,'_dd_mm_yyyy_HH_MM_SS')];
    saveas(H,[SAVE_NAME,'.fig'])
end

end


%Code to accompany the paper:
%"Modelling persistence of motion in a crowded environment: the diffusive
%limit of excluding velocity-jump processes"
%by Enrico Gavagnin and Christian A. Yates

%Created 06/10/2017
%By Enrico Gavagnin
%email: e.gavagnin@bath.ac.uk
%%
% This function is designed to simulate the agent-based model with
% parameters received as input and it produces as output the matrix
%
function [Rx_mean,Lx_mean,Ux_mean,Dx_mean]=ABM(type,x_size,y_size,IC,P_m,phi,P_r,v,M,T_final,d)


%% INPUT
% type: type of agent interaction to implement
% x_size: number of column of the lattice
% y_sites: number of rows of the lattice
% IC: inital condition parameter
% P_m: motility rate
% phi: persistence parameter
% P_r: reorienting rate
% v: velocity
% M: number of rrealisations of the simulation
% T_final: final time for each simulation
% d: inital density

%% OUTPUT
% Rx/Lx/Ux/Dx_mean: matrix counting the number of Right/Left/Up/Down-type
% particles average over many realisations

%Initiate the matrices of the densities
Rx=zeros(y_size,x_size);
Lx=zeros(y_size,x_size);
Ux=zeros(y_size,x_size);
Dx=zeros(y_size,x_size);

% Repeats loop
for m=1:M
    %Initiate time variable
    t=0;
    
    %Initiate the matrix of the lattice
    L=Initial_Conditions(IC,x_size,y_size,d);
    
    %Compute the total number of agents
    tot_agts=nnz(L);
    
    %Initiate the matrix of the coordinates of all the agents
    C=zeros(tot_agts,3);
    [C(:,1),C(:,2),C(:,3)]=find(L);
    
    
    
    %% Time loop
    while t<T_final
        
        %Generate the time of the next event
        tau=log(1/rand)/((P_m+P_r)*tot_agts);
        
        %Choose an agent to which an event applies
        agt=ceil(tot_agts*rand);
        
        % Polarisations legend:
        % j=    2: Right
        %      -2: Left
        %     0.5: Up
        %    -0.5: Down
        
        %Initialise the vector with the possible polarisations (according
        %to the dimension of the lattice)
        
        if y_size~=1
            pol_vector=[2,-2,0.5,-0.5];
        else
            pol_vector=[2,-2];
        end
        
        %Decide the event that occurs between reorient and movement
        %% Reorient
        if rand<P_r/(P_r+P_m)
            
            %Uniformly choose the new specie for the particle
            C(agt,3)=pol_vector(ceil(size(pol_vector,2)*rand));
            L(C(agt,1),C(agt,2))=C(agt,3);
            %% Movement
        else
            
            % Choose the direction of movement according to the
            % polarisation  of the selected agent
            
            rn=rand;
            % If the simulation is in two dimensions
            if y_size~=1
                
                if rn<(1+phi)/4
                    j=C(agt,3);
                elseif rn<0.5
                    j=-C(agt,3);
                elseif rn<3/4
                    j=1/C(agt,3);
                else
                    j=-1/C(agt,3);
                end
                %if the simulation is in one dimension
            else
                if rn<(1+phi)/2
                    j=C(agt,3);
                else
                    j=-C(agt,3);
                end
            end
            
            
            % Check if the movement is allowed, accordign to the type of
            % interaction selected and that case update the agent's position
            
            
            
            if type==0 %TYPE 0
                %%
                % If the event is RIGHT
                if (j==2)&&(C(agt,2)<x_size-v)
                    C(agt,2)=C(agt,2)+v;
                end
                
                % If the event is LEFT
                if (j==-2)&&(C(agt,2)>v)
                    C(agt,2)=C(agt,2)-v;
                end
                
                % If the event is UP
                if (j==0.5)
                    if (C(agt,1)>v)
                        C(agt,1)=C(agt,1)-v;
                    else
                        C(agt,1)=y_size-(v-C(agt,1));
                    end
                end
                
                % If the event is DOWN
                if (j==-0.5)
                    if (C(agt,1)<y_size-v)
                        C(agt,1)=C(agt,1)+v;
                    else
                        C(agt,1)=C(agt,1)+v-y_size+1;
                    end
                end
                
            elseif type==1 %TYPE 1
                %%
                % If the event is RIGHT
                if (j==2)&&(C(agt,2)<x_size-v)&&(L(C(agt,1),C(agt,2)+v)==0)
                    L(C(agt,1),C(agt,2)+v)=C(agt,3);
                    L(C(agt,1),C(agt,2))=0;
                    C(agt,2)=C(agt,2)+v;
                end%if
                
                % If the event is LEFT
                if (j==-2)&&(C(agt,2)>v)&&(L(C(agt,1),C(agt,2)-v)==0)
                    L(C(agt,1),C(agt,2)-v)=C(agt,3);
                    L(C(agt,1),C(agt,2))=0;
                    C(agt,2)=C(agt,2)-v;
                end %if
                
                % If the event is UP
                if (j==0.5)
                    if (C(agt,1)>v)
                        if (L(C(agt,1)-v,C(agt,2))==0)
                            L(C(agt,1)-v,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)-v;
                        end
                    else
                        if (L(C(agt,1)-v+y_size,C(agt,2))==0)
                            L(C(agt,1)+y_size-v,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)+y_size-v;
                        end
                    end
                end %if
                
                % If the event is DOWN
                if (j==-0.5)
                    if C(agt,1)<=y_size-v
                        if (L(C(agt,1)+v,C(agt,2))==0)
                            L(C(agt,1)+v,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)+v;
                        end
                    else
                        if (L(C(agt,1)+v-y_size,C(agt,2))==0)
                            L(C(agt,1)+v-y_size,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)+v-y_size;
                        end
                    end
                end %if
                
                
            elseif type==2 %TYPE 2
                %%
                % If the event is RIGHT
                if (j==2)&&(C(agt,2)<x_size-v)
                    if (nnz(L(C(agt,1),C(agt,2)+1:C(agt,2)+v))==0)
                        L(C(agt,1),C(agt,2)+v)=C(agt,3);
                        L(C(agt,1),C(agt,2))=0;
                        C(agt,2)=C(agt,2)+v;
                    end
                end%if
                
                % If the event is LEFT
                if (j==-2)&&(C(agt,2)>v)
                    
                    if (nnz(L(C(agt,1),C(agt,2)-v:C(agt,2)-1))==0)
                        L(C(agt,1),C(agt,2)-v)=C(agt,3);
                        L(C(agt,1),C(agt,2))=0;
                        C(agt,2)=C(agt,2)-v;
                    end
                end %if
                
                % If the event is UP
                if (j==0.5)
                    if (C(agt,1)>v)
                        if (nnz(L(C(agt,1)-v:C(agt,1)-1,C(agt,2)))==0)
                            L(C(agt,1)-v,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)-v;
                        end
                    else
                        if (nnz([L(1:C(agt,1)-1,C(agt,2));L(end+C(agt,1)-v:end,C(agt,2))])==0)
                            L(y_size+C(agt,1)-v,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=y_size+C(agt,1)-v;
                        end
                        
                    end
                end %if
                
                % If the event is DOWN
                if (j==-0.5)
                    if C(agt,1)<=y_size-v
                        if (nnz(L(C(agt,1)+1:C(agt,1)+v,C(agt,2)))==0)
                            L(C(agt,1)+v,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)+v;
                        end
                    else
                        if (nnz([L(C(agt,1)+1:end,C(agt,2));L(1:C(agt,1)+v-y_size,C(agt,2))])==0)
                            L(C(agt,1)+v-y_size,C(agt,2))=C(agt,3);
                            L(C(agt,1),C(agt,2))=0;
                            C(agt,1)=C(agt,1)+v-y_size;
                        end
                    end
                    
                end %if
                
            elseif type==3 %TYPE 3
                %%
                scout=0; %inspector index
                insp_value=0;
                
                % If the event is RIGHT
                if (j==2)
                    av_steps=x_size-C(agt,2); %Boundary conditions
                    
                    if av_steps~=0
                        %Find the farest available site up to N
                        while (av_steps>0)&&(scout<v)&&(insp_value==0)
                            scout=scout+1;
                            av_steps=av_steps-1;
                            insp_value=L(C(agt,1),C(agt,2)+scout);
                        end% while
                        
                        %Check if the while loop found an available site
                        if (insp_value==0) || (scout==v)
                            scout=scout-((insp_value~=0)*(scout==v));
                            
                            %Fix the inspector index if it is N and occupied
                            scout=scout-((insp_value~=0)*(scout==v));
                            
                            %Update lattice matrix
                            L(C(agt,1),C(agt,2))=0;
                            L(C(agt,1),C(agt,2)+scout)=C(agt,3);
                            
                            
                            %Update cooordinates matrix
                            C(agt,2)=C(agt,2)+scout;
                        end%if
                    end%if
                end%if
                
                % If the event is LEFT
                if (j==-2)
                    av_steps=C(agt,2)-1; %Boundary conditions
                    
                    if av_steps~=0
                        %Find the farest available site up to N
                        while (av_steps>0)&&(scout<v)&&(insp_value==0)
                            scout=scout+1;
                            av_steps=av_steps-1;
                            insp_value=L(C(agt,1),C(agt,2)-scout);
                        end% while
                        
                        %Check if the while loop found an available site
                        if (insp_value==0) || (scout==v)
                            
                            %Fix the inspector index if it is N and occupied
                            scout=scout-((insp_value~=0)*(scout==v));
                            
                            %Update lattice matrix
                            L(C(agt,1),C(agt,2))=0;
                            
                            L(C(agt,1),C(agt,2)-scout)=C(agt,3);
                            
                            %Update cooordinates matrix
                            C(agt,2)=C(agt,2)-scout;
                        end%if
                    end%if
                end %if
                
                % If the event is UP
                if (j==0.5)
                    
                    
                    %Find the farest available site up to N
                    while (scout<v)&&(insp_value==0)
                        scout=scout+1;
                        insp_value=L(mod(C(agt,1)-scout-1,y_size)+1,C(agt,2));
                    end% while
                    
                    %Check if the while loop found an available site
                    if (insp_value==0) || (scout==v)
                        %Fix the inspector index if it is N and occupied
                        scout=scout-((insp_value~=0)*(scout==v));
                        
                        if C(agt,1)-scout>0
                            %Update lattice matrix
                            L(C(agt,1),C(agt,2))=0;
                            L(C(agt,1)-scout,C(agt,2))=C(agt,3);
                            
                            %Update cooordinates matrix
                            C(agt,1)=C(agt,1)-scout;
                        else
                            %Update lattice matrix
                            L(C(agt,1),C(agt,2))=0;
                            L(C(agt,1)-scout+y_size,C(agt,2))=C(agt,3);
                            
                            %Update cooordinates matrix
                            C(agt,1)=C(agt,1)-scout+y_size;
                        end
                        
                    end%if
                end %if
                
                % If the event is DOWN
                if (j==-0.5)
                    %Find the farest available site up to N
                    while (scout<v)&&(insp_value==0)
                        scout=scout+1;
                        insp_value=L(mod(C(agt,1)+scout-1,y_size)+1,C(agt,2));
                    end% while
                    
                    %Check if the while loop found an available site
                    if (insp_value==0) || (scout==v)
                        %Fix the inspector index if it is N and occupied
                        scout=scout-((insp_value~=0)*(scout==v));
                        
                        if C(agt,1)+scout<=y_size
                            %Update lattice matrix
                            L(C(agt,1),C(agt,2))=0;
                            L(C(agt,1)+scout,C(agt,2))=C(agt,3);
                            
                            %Update cooordinates matrix
                            C(agt,1)=C(agt,1)+scout;
                        else
                            %Update lattice matrix
                            L(C(agt,1),C(agt,2))=0;
                            L(C(agt,1)+scout-y_size,C(agt,2))=C(agt,3);
                            
                            %Update cooordinates matrix
                            C(agt,1)=C(agt,1)+scout-y_size;
                        end
                    end%if
                    
                end
            end
            
            %% Update time
            t=t+tau;
        end
        
    end %for
    %Update sum of relisations
        if type==0
            R_idx=C(:,3)==2;
            Rx=Rx+sparse(C(R_idx,1),C(R_idx,2),ones(nnz(R_idx),1),y_size,x_size);
            L_idx=C(:,3)==-2;
            Lx=Lx+sparse(C(L_idx,1),C(L_idx,2),ones(nnz(L_idx),1),y_size,x_size);
            U_idx=C(:,3)==0.5;
            Ux=Ux+sparse(C(U_idx,1),C(U_idx,2),ones(nnz(U_idx),1),y_size,x_size);
            D_idx=C(:,3)==-0.5;
            Dx=Dx+sparse(C(D_idx,1),C(D_idx,2),ones(nnz(D_idx),1),y_size,x_size);
        else
            Rx=Rx+(L==2);
            Lx=Lx+(L==-2);
            Ux=Ux+(L==0.5);
            Dx=Dx+(L==-0.5);
        end
end
% Normalise the number of agts per site by the number of realisations
Rx_mean=Rx/(M);
Lx_mean=Lx/(M);
Ux_mean=Ux/(M);
Dx_mean=Dx/(M);
end
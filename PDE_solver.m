%Code to accompany the paper:
%"Modelling persistence of motion in a crowded environment: the diffusive
%limit of excluding velocity-jump processes"
%by Enrico Gavagnin and Christian A. Yates

%Created 06/10/2017
%By Enrico Gavagnin
%email: e.gavagnin@bath.ac.uk
%%
% This function is designed to solve the column-averaged PDEs using
% implicit Euler method with Picard iterations.

function [R_n, L_n, U_n, D_n]=PDE_solver(type,x_size,u_IC,T_final,P_m,phi,P_r,v)


%% INPUT
% type: type of agent interaction to implement
% x_size: number of column of the lattice
% u_IC: inital condition for the numerical solution
% T_final: final time for each simulation
% P_m: motility rate
% phi: persistence parameter
% P_r: reorienting rate
% v: velocity


%% OUTPUT
% Rn/Ln/Un/Dn: vectors of the numerical solution of Right/Left/Up/Down
% column-averaged density at time T_final



% PICARD ITERATION
%Given the that the discretisation of the PDE is written in the form:
%                       A(u)*u=b(u)
%At each iteration we solve the linear system
%                       A(u-)*u+=b(u-)
%and we then update u-=u+ iteratively until we have dist(u-,u+)< tolerance.
%The process starts with u- equal to the most recent numerical solution
%computed u, and u+ equal to the vecor all of all ones to make sure we
%enter in the while loop at the first step.
%
%We define A as a sparse matrix.

%% Initialisation

%Define spacial step for the numerical solution
delta_x=1*10^(0);

%Define time step for the numerical solution
delta_t=2*10^(-1);

%Define the tolerance parameter for the Picard method
toll=10^(-3);

%Define maximum number of iterations for Picard method
max_iter=100;

%Initiate time extremal value
T_final_PDE=T_final/delta_t;

%Initiate numerical solution
%u=[ u_R, u_L, u_U, ,u_D]
u=u_IC;

%Store the spacial size of the numerical solution
row=x_size/delta_x;

%Preallocate vector coord which store the non zeros elements of the
%matrix A for the Picard Iteration.
coord=[];

if type==0
    
    %Define three constants diff_const, adv_const and reac_const
    diff_const=P_m*v^2*delta_t/(4*delta_x^2);
    adv_const=P_m*v*phi*delta_t/(4*delta_x);
    react_const=delta_t*P_r/4;
    
    
    %% Time loop
    for k=1:T_final_PDE
        
        
        % Initialisation variables for each Picard Iteration
        u_plus=ones(size(u));
        iter=1;
        u_minus=u;
        u_temp=u;
        
        %% Iteration loop
        while iter<max_iter && norm(u_plus-u_minus)>toll
            
            % Update the vector u-
            u_minus=u_temp;
            
            %Define a variable counter that keeps track of the element stored
            %in coord
            counter=0;
            
            %% Update the matrix A
            
            %First quarter of the matrix A
            for i=1:row
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    %Update the vector of the non-xeros elemens of A
                    coord(counter+1:counter+6,:)=...
                        [i, i, 1+2*diff_const+3*react_const;     %Coefficient terms of R_{i,i}
                        i, i-1, -diff_const-adv_const;           %Coefficient terms of R_{i,i-1}
                        i, i+1, -diff_const+adv_const;           %Coefficient terms of R_{i,i+1}
                        i, i+row, -react_const;                  %Coefficient terms of L_{i,j}
                        i, i+2*row, -react_const;                %Coefficient terms of U_{i,j}
                        i, i+3*row, -react_const;                %Coefficient terms of D_{i,j}
                        ];
                    counter=counter+6;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            %Second quarter of the matrix A
            for i=row+1:2*row
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    %Update the vector of the non-xeros elemens of A
                    coord(counter+1:counter+6,:)=...
                        [i, i, 1+2*diff_const+3*react_const;     %Coefficient terms of L_{i,i}
                        i, i-1, -diff_const+adv_const;           %Coefficient terms of L_{i,i-1}
                        i, i+1, -diff_const-adv_const;           %Coefficient terms of L_{i,i+1}
                        i, i-row, -react_const;                  %Coefficient terms of R_{i,j}
                        i, i+row, -react_const;                  %Coefficient terms of U_{i,j}
                        i, i+2*row, -react_const;                %Coefficient terms of D_{i,j}
                        ];
                    counter=counter+6;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            %Third quarter of the matrix A
            for i=2*row+1:3*row
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    %Update the vector of the non-xeros elemens of A
                    coord(counter+1:counter+6,:)=...
                        [i, i, 1+2*diff_const+3*react_const;       %Coefficient terms of U_{i,i}
                        i, i-1, -diff_const;                       %Coefficient terms of U_{i,i-1}
                        i, i+1, -diff_const;                       %Coefficient terms of U_{i+1,i}
                        i, i-2*row, -react_const;                  %Coefficient terms of R_{i,j}
                        i, i-row, -react_const;                    %Coefficient terms of L_{i,j}
                        i, i+row, -react_const;                    %Coefficient terms of D_{i,j}
                        ];
                    counter=counter+6;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            %Third quarter of the matrix A
            for i=3*row+1:4*row
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    coord(counter+1:counter+6,:)=...
                        [i, i, 1+2*diff_const+3*react_const;       %Coefficient terms of D_{i,i}
                        i, i-1, -diff_const;                       %Coefficient terms of D_{i,i-1}
                        i, i+1, -diff_const;                       %Coefficient terms of D_{i+1,i}
                        i, i-3*row, -react_const;                  %Coefficient terms of R_{i,j}
                        i, i-2*row, -react_const;                  %Coefficient terms of L_{i,j}
                        i, i-row, -react_const;                   %Coefficient terms of U_{i,j}
                        ];
                    counter=counter+6;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            
            %% Solution of the linear system
            
            % Split the solution for the different polarisations
            R=u(1:row)';
            L=u(row+1:2*row)';
            U=u(2*row+1:3*row)';
            D=u(3*row+1:4*row)';
            
            %Define the vector b of the linear system
            %We force certain elements of the vector b to impose the boundary condition
            b=[0;R(2:end-1);0;0;L(2:end-1);0;0;U(2:end-1);0;0;D(2:end-1);0]; 
            
            %Build the sparse matrix A
            A=sparse(coord(:,1),coord(:,2),coord(:,3),size(u_IC,2),size(u_IC,2));
            
            %Solve the linear system
            u_plus=A\b;
            u_plus=u_plus';
            %Update auxiliary vector
            u_temp=u_plus;
            
            %Update number of iteration
            iter=iter+1;
            
        end %while
        %%
        
        %Display an error message if the tolerance has not been reached
        if iter==max_iter
            error('Picard Iteration reaches maximum iteration number');
        end % if
        
        % Update the Backward Euler time iteration
        u=u_plus;
        
    end % for (TIME LOOP)
    
elseif type==1
    
    %Define three constants diff_const, adv_const and reac_const
    diff_const=P_m*v^2*delta_t/(4*delta_x^2);
    adv_const=P_m*v*phi*delta_t/(4*delta_x);
    react_const=delta_t*P_r/4;
    
    %% Time loop
    for k=1:T_final_PDE
        
        % Initialisation variables for each Picard Iteration
        u_plus=ones(size(u));
        iter=1;
        u_minus=u;
        u_temp=u;
        
        %% Iteration loop
        while iter<max_iter && norm(u_plus-u_minus)>toll
            % Update the vector u-
            u_minus=u_temp;
            
            
            %Split the u vector in the two corresponding halfs
            R=u_minus(1:row);
            L=u_minus(row+1:2*row);
            U=u_minus(2*row+1:3*row);
            D=u_minus(3*row+1:4*row);
            C=R+L+U+D;
            
            
            %Define a variable counter that keeps track of the element stored
            %in coord
            
            counter=0;
            
            % fisrt quorter of A rows
            for i=1:row
                if mod(i,row)>1
                    %Indeces corresponding to non-boundary sites
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+2*diff_const*(R(i)+1-C(i))+3*react_const;                %Coefficient terms of R_{i}
                        i, i-1, -diff_const*(R(i)+1-C(i))-adv_const*(1-C(i-1));           %Coefficient terms of R_{i-1}
                        i, i+1, -diff_const*(R(i)+1-C(i))+adv_const*(1-C(i+1));           %Coefficient terms of R_{i+1}
                        
                        i, i+row-1, -diff_const*R(i);                                     %Coefficient terms of L_{i-1}
                        i, i+row, 2*diff_const*R(i)-react_const;                          %Coefficient terms of L_{i}
                        i, i+row+1, -diff_const*R(i);                                     %Coefficient terms of L_{i+1}
                        
                        i, i+2*row-1, -diff_const*R(i);                                   %Coefficient terms of U_{i-1}
                        i, i+2*row, 2*diff_const*R(i)-react_const;                        %Coefficient terms of U_{i,j}
                        i, i+2*row+1, -diff_const*R(i);                                   %Coefficient terms of U_{i+1}
                        
                        i, i+3*row-1, -diff_const*R(i);                                   %Coefficient terms of D_{i-1}
                        i, i+3*row, 2*diff_const*R(i)-react_const;                        %Coefficient terms of D_{i,j}
                        i, i+3*row+1, -diff_const*R(i);                                   %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            % second quorter of A rows
            for i=row+1:2*row
                ii=i-row;
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+2*diff_const*(L(ii)+1-C(ii))+3*react_const;                 %Coefficient terms of L_{i}
                        i, i-1, -diff_const*(L(ii)+1-C(ii))+adv_const*(1-C(ii-1));           %Coefficient terms of L_{i-1}
                        i, i+1, -diff_const*(L(ii)+1-C(ii))-adv_const*(1-C(ii+1));           %Coefficient terms of L_{i+1}
                        
                        i, i-row-1, -diff_const*L(ii);                                       %Coefficient terms of R_{i-1}
                        i, i-row, 2*diff_const*L(ii)-react_const;                            %Coefficient terms of R_{i}
                        i, i-row+1, -diff_const*L(ii);                                       %Coefficient terms of R_{i+1}
                        
                        i, i+row-1, -diff_const*L(ii);                                       %Coefficient terms of U_{i-1}
                        i, i+row, 2*diff_const*L(ii)-react_const;                            %Coefficient terms of U_{i,j}
                        i, i+row+1, -diff_const*L(ii);                                       %Coefficient terms of U_{i+1}
                        
                        i, i+2*row-1, -diff_const*L(ii);                                     %Coefficient terms of D_{i-1}
                        i, i+2*row, 2*diff_const*L(ii)-react_const;                          %Coefficient terms of D_{i,j}
                        i, i+2*row+1, -diff_const*L(ii);                                     %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            % Third quorter of A rows
            for i=2*row+1:3*row
                ii=i-2*row;
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+2*diff_const*(U(ii)+1-C(ii))+3*react_const;       %Coefficient terms of U_{i}
                        i, i-1, -diff_const*(U(ii)+1-C(ii));                       %Coefficient terms of U_{i-1}
                        i, i+1, -diff_const*(U(ii)+1-C(ii));                       %Coefficient terms of U_{i+1}
                        
                        i, i-2*row-1, -diff_const*U(ii);                           %Coefficient terms of R_{i-1}
                        i, i-2*row, 2*diff_const*U(ii)-react_const;                %Coefficient terms of R_{i}
                        i, i-2*row+1, -diff_const*U(ii);                           %Coefficient terms of R_{i+1}
                        
                        i, i-row-1, -diff_const*U(ii);                             %Coefficient terms of L_{i-1}
                        i, i-row, 2*diff_const*U(ii)-react_const;                  %Coefficient terms of L_{i}
                        i, i-row+1, -diff_const*U(ii);                             %Coefficient terms of L_{i+1}
                        
                        i, i+row-1, -diff_const*U(ii);                             %Coefficient terms of D_{i-1}
                        i, i+row, 2*diff_const*U(ii)-react_const;                  %Coefficient terms of D_{i}
                        i, i+row+1, -diff_const*U(ii);                             %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                else
                    
                    %Indeces corresponding to boundary sites
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            % Fourth quorter of A rows
            for i=3*row+1:4*row
                ii=i-3*row;
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+2*diff_const*(D(ii)+1-C(ii))+3*react_const;       %Coefficient terms of D_{i}
                        i, i-1, -diff_const*(D(ii)+1-C(ii));                       %Coefficient terms of D_{i-1}
                        i, i+1, -diff_const*(D(ii)+1-C(ii));                       %Coefficient terms of D_{i+1}
                        
                        i, i-3*row-1, -diff_const*D(ii);                           %Coefficient terms of R_{i-1}
                        i, i-3*row, 2*diff_const*D(ii)-react_const;                %Coefficient terms of R_{i}
                        i, i-3*row+1, -diff_const*D(ii);                           %Coefficient terms of R_{i+1}
                        
                        i, i-2*row-1, -diff_const*D(ii);                           %Coefficient terms of L_{i-1}
                        i, i-2*row, 2*diff_const*D(ii)-react_const;                %Coefficient terms of L_{i}
                        i, i-2*row+1, -diff_const*D(ii);                           %Coefficient terms of L_{i+1}
                        
                        i, i-row-1, -diff_const*D(ii);                             %Coefficient terms of U_{i-1}
                        i, i-row, 2*diff_const*D(ii)-react_const;                  %Coefficient terms of U_{i}
                        i, i-row+1, -diff_const*D(ii);                             %Coefficient terms of U_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            %% Solution of the linear system
            
            % Split the solution for the different polarisations
            R=u(1:row)';
            L=u(row+1:2*row)';
            U=u(2*row+1:3*row)';
            D=u(3*row+1:4*row)';
            
            %Define the vector b of the linear system
            %We force certain elements of the vector b to impose the boundary condition
            b=[0;R(2:end-1);0;0;L(2:end-1);0;0;U(2:end-1);0;0;D(2:end-1);0]; 
            
            %Build the sparse matrix A
            A=sparse(coord(:,1),coord(:,2),coord(:,3),size(u_IC,2),size(u_IC,2));
            
            %Solve the linear system
            u_plus=A\b;
            u_plus=u_plus';
            %Update auxiliary vector
            u_temp=u_plus;
            
            %Update number of iteration
            iter=iter+1;
            
        end %while
        %%
        
        %Display an error message if the tolerance has not been reached
        if iter==max_iter
            error('Picard Iteration reaches maximum iteration number');
        end % if
        
        % Update the Backward Euler time iteration
        u=u_plus;
        
    end % for (TIME LOOP)
    
elseif type==2
    
    %Define three constants diff_const, adv_const and reac_const
    diff_const=P_m*v^2*delta_t/(4*delta_x^2);
    adv_const=P_m*v*phi*delta_t/(4*delta_x);
    react_const=delta_t*P_r/4;
    
    
    %% Time loop
    for k=1:T_final_PDE
        
        % Initialisation variables for each Picard Iteration
        u_plus=ones(size(u));
        iter=1;
        u_minus=u;
        u_temp=u;
        
        %% Iteration loop
        while iter<max_iter && norm(u_plus-u_minus)>toll
            % Update the vector u-
            u_minus=u_temp;
            
            
            %Split the u vector in the two corresponding halfs
            R=u_minus(1:row);
            L=u_minus(row+1:2*row);
            U=u_minus(2*row+1:3*row);
            D=u_minus(3*row+1:4*row);
            C=R+L+U+D;
            
            
            %Define a variable counter that keeps track of the element stored
            %in coord
            counter=0;
            
            %Loop trough the rows
            
            %First quorter
            for i=1:row
                %Define the diffusion coefficient
                diff_coeff=diff_const*(1-C(i))^(v-2);
                
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    
                    C_diff=(C(i+1)-C(i-1))/2;
                    
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_coeff*(2*(1-C(i))^2+2*R(i)*(1-C(i))+(v-1)*C_diff^2)+3*react_const;                %Coefficient terms of R_{i}
                        i, i-1, diff_coeff*(-(v-1)*(1-C(i))/2*C_diff-(1-C(i))^2-R(i)*(1-C(i)))-adv_const*(1-C(i-1))^v;  %Coefficient terms of R_{i-1}
                        i, i+1, diff_coeff*((v-1)*(1-C(i))/2*C_diff-(1-C(i))^2-R(i)*(1-C(i)))+adv_const*(1-C(i+1))^v;   %Coefficient terms of R_{i+1}
                        
                        i, i+row-1, -diff_coeff*R(i)*(1-C(i));                   %Coefficient terms of L_{i-1}
                        i, i+row, 2*diff_coeff*R(i)*(1-C(i))-react_const;        %Coefficient terms of L_{i}
                        i, i+row+1, -diff_coeff*R(i)*(1-C(i));                   %Coefficient terms of L_{i+1}
                        
                        i, i+2*row-1, -diff_coeff*R(i)*(1-C(i));                  %Coefficient terms of U_{i-1}
                        i, i+2*row, 2*diff_coeff*R(i)*(1-C(i))-react_const;       %Coefficient terms of U_{i,j}
                        i, i+2*row+1, -diff_coeff*R(i)*(1-C(i));                  %Coefficient terms of U_{i+1}
                        
                        i, i+3*row-1, -diff_coeff*R(i)*(1-C(i));                  %Coefficient terms of D_{i-1}
                        i, i+3*row, 2*diff_coeff*R(i)*(1-C(i))-react_const;       %Coefficient terms of D_{i,j}
                        i, i+3*row+1, -diff_coeff*R(i)*(1-C(i));                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            % Second quorter
            for i=row+1:2*row
                
                ii=i-row;
                % Define the diffusion coefficient
                diff_coeff=diff_const*(1-C(ii))^(v-2);
                
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    C_diff=(C(ii+1)-C(ii-1))/2;
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_coeff*(2*(1-C(ii))^2+2*L(ii)*(1-C(ii))+(v-1)*C_diff^2)+3*react_const;                      %Coefficient terms of L_{i}
                        i, i-1, diff_coeff*(-(v-1)*(1-C(ii))/2*C_diff-(1-C(ii))^2-L(ii)*(1-C(ii)))+adv_const*(1-C(ii-1))^v;      %Coefficient terms of L_{i-1}
                        i, i+1, diff_coeff*((v-1)*(1-C(ii))/2*C_diff-(1-C(ii))^2-L(ii)*(1-C(ii)))-adv_const*(1-C(ii+1))^v;       %Coefficient terms of L_{i+1}
                        
                        i, i-row-1, -diff_coeff*L(ii)*(1-C(ii));                   %Coefficient terms of R_{i-1}
                        i, i-row, 2*diff_coeff*L(ii)*(1-C(ii))-react_const;        %Coefficient terms of R_{i}
                        i, i-row+1, -diff_coeff*L(ii)*(1-C(ii));                   %Coefficient terms of R_{i+1}
                        
                        i, i+row-1, -diff_coeff*L(ii)*(1-C(ii));                  %Coefficient terms of U_{i-1}
                        i, i+row, 2*diff_coeff*L(ii)*(1-C(ii))-react_const;       %Coefficient terms of U_{i,j}
                        i, i+row+1, -diff_coeff*L(ii)*(1-C(ii));                  %Coefficient terms of U_{i+1}
                        
                        i, i+2*row-1, -diff_coeff*L(ii)*(1-C(ii));                  %Coefficient terms of D_{i-1}
                        i, i+2*row, 2*diff_coeff*L(ii)*(1-C(ii))-react_const;       %Coefficient terms of D_{i,j}
                        i, i+2*row+1, -diff_coeff*L(ii)*(1-C(ii));                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            % Third quorter
            for i=2*row+1:3*row
                ii=i-2*row;
                %Define the diffusion coefficient
                diff_coeff=diff_const*(1-C(ii))^(v-2);
                
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    C_diff=(C(ii+1)-C(ii-1))/2;
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_coeff*(2*(1-C(ii))^2+2*U(ii)*(1-C(ii))+(v-1)*C_diff^2)+3*react_const;       %Coefficient terms of U_{i}
                        i, i-1, diff_coeff*(-(v-1)*(1-C(ii))/2*C_diff-(1-C(ii))^2-U(ii)*(1-C(ii)));               %Coefficient terms of U_{i-1}
                        i, i+1, diff_coeff*((v-1)*(1-C(ii))/2*C_diff-(1-C(ii))^2-U(ii)*(1-C(ii)));                %Coefficient terms of U_{i+1}
                        
                        i, i-2*row-1, -diff_coeff*U(ii)*(1-C(ii));                   %Coefficient terms of R_{i-1}
                        i, i-2*row, 2*diff_coeff*U(ii)*(1-C(ii))-react_const;        %Coefficient terms of R_{i}
                        i, i-2*row+1, -diff_coeff*U(ii)*(1-C(ii));                   %Coefficient terms of R_{i+1}
                        
                        i, i-row-1, -diff_coeff*U(ii)*(1-C(ii));                  %Coefficient terms of L_{i-1}
                        i, i-row, 2*diff_coeff*U(ii)*(1-C(ii))-react_const;       %Coefficient terms of L_{i}
                        i, i-row+1, -diff_coeff*U(ii)*(1-C(ii));                  %Coefficient terms of L_{i+1}
                        
                        i, i+row-1, -diff_coeff*U(ii)*(1-C(ii));                  %Coefficient terms of D_{i-1}
                        i, i+row, 2*diff_coeff*U(ii)*(1-C(ii))-react_const;       %Coefficient terms of D_{i}
                        i, i+row+1, -diff_coeff*U(ii)*(1-C(ii));                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            % Fourth quorter
            for i=3*row+1:4*row
                ii=i-3*row;
                %Define the diffusion coefficient
                diff_coeff=diff_const*(1-C(ii))^(v-2);
                
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    C_diff=(C(ii+1)-C(ii-1))/2;
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_coeff*(2*(1-C(ii))^2+2*D(ii)*(1-C(ii))+(v-1)*C_diff^2)+3*react_const;       %Coefficient terms of D_{i}
                        i, i-1, diff_coeff*(-(v-1)*(1-C(ii))/2*C_diff-(1-C(ii))^2-D(ii)*(1-C(ii)));               %Coefficient terms of D_{i-1}
                        i, i+1,diff_coeff*((v-1)*(1-C(ii))/2*C_diff-(1-C(ii))^2-D(ii)*(1-C(ii)));                 %Coefficient terms of D_{i+1}
                        
                        i, i-3*row-1, -diff_coeff*D(ii)*(1-C(ii));                   %Coefficient terms of R_{i-1}
                        i, i-3*row, 2*diff_coeff*D(ii)*(1-C(ii))-react_const;        %Coefficient terms of R_{i}
                        i, i-3*row+1, -diff_coeff*D(ii)*(1-C(ii));                   %Coefficient terms of R_{i+1}
                        
                        i, i-2*row-1, -diff_coeff*D(ii)*(1-C(ii));                  %Coefficient terms of L_{i-1}
                        i, i-2*row, 2*diff_coeff*D(ii)*(1-C(ii))-react_const;       %Coefficient terms of L_{i}
                        i, i-2*row+1, -diff_coeff*D(ii)*(1-C(ii));                  %Coefficient terms of L_{i+1}
                        
                        i, i-row-1, -diff_coeff*D(ii)*(1-C(ii));                  %Coefficient terms of U_{i-1}
                        i, i-row, 2*diff_coeff*D(ii)*(1-C(ii))-react_const;       %Coefficient terms of U_{i}
                        i, i-row+1, -diff_coeff*D(ii)*(1-C(ii));                  %Coefficient terms of U_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            %% Solution of the linear system
            
            % Split the solution for the different polarisations
            R=u(1:row)';
            L=u(row+1:2*row)';
            U=u(2*row+1:3*row)';
            D=u(3*row+1:4*row)';
            
            %Define the vector b of the linear system
            %We force certain elements of the vector b to impose the boundary condition
            b=[0;R(2:end-1);0;0;L(2:end-1);0;0;U(2:end-1);0;0;D(2:end-1);0]; 
            
            %Build the sparse matrix A
            A=sparse(coord(:,1),coord(:,2),coord(:,3),size(u_IC,2),size(u_IC,2));
            
            %Solve the linear system
            u_plus=A\b;
            u_plus=u_plus';
            %Update auxiliary vector
            u_temp=u_plus;
            
            %Update number of iteration
            iter=iter+1;
            
        end %while
        %%
        
        %Display an error message if the tolerance has not been reached
        if iter==max_iter
            error('Picard Iteration reaches maximum iteration number');
        end % if
        
        % Update the Backward Euler time iteration
        u=u_plus;
        
    end % for (TIME LOOP)
    
elseif type==3
    
    %Define three constants diff_const, adv_const and reac_const
    diff_const=P_m*delta_t/(4*delta_x^2);
    adv_const=P_m*phi*delta_t/(4*delta_x);
    react_const=delta_t*P_r/4;
    
    if v>3
        error('the numerical solution for the type 3 interaction is implemented only for v<4')
    end
    v_2=0;
    v_3=0;
    if v>1
        v_2=1;
    end
    if v>2
        v_3=1;
    end
    
    for k=1:T_final_PDE
        
        % Initialisation variables for each Picard Iteration
        u_plus=ones(size(u));
        iter=1;
        u_minus=u;
        u_temp=u;
        
        %% Iteration loop
        while iter<max_iter && norm(u_plus-u_minus)>toll
            
            % Update the vector u-
            u_minus=u_temp;
            
            
            %Split the u vector in the two corresponding halfs
            R=u_minus(1:row);
            L=u_minus(row+1:2*row);
            U=u_minus(2*row+1:3*row);
            D=u_minus(3*row+1:4*row);
            C=R+L+U+D;
            
            
            %Define a variable counter that keeps track of the element stored
            %in coord
            
            counter=0;
            
            %Loop trough the rows
            
            % First quorter
            for i=1:row 
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    
                    %Define the finite difference of the total density
                    %C_diff
                    C_diff=(C(i+1)-C(i-1))/2;
                    if C(i+1)~=0
                        C_adv_plus=(1-C(i+1))*((1-C(i+1))^v-1)/(C(i+1));
                    else
                        C_adv_plus=0;
                    end
                    if C(i-1)~=0
                        C_adv_minus=(1-C(i-1))*((1-C(i-1))^v-1)/(C(i-1));
                    else
                        C_adv_minus=0;
                    end
                    
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_const*(2*(1-C(i))+2*R(i)+v_2*6*(1-C(i))^2+v_3*(10*(1-C(i))^3-6*(1-C(i))*C_diff^2+3*(1-C(i))^2*(C(i+1)-2*C(i)+C(i-1))))+3*react_const;       %Coefficient terms of R_{i}
                        i, i-1, diff_const*(-(1-C(i))-R(i)+v_2*(-3*(1-C(i))*C_diff-3*(1-C(i))^2)+v_3*(-9*(1-C(i))^2*C_diff-5*(1-C(i))^3))+adv_const*C_adv_minus;                  %Coefficient terms of R_{i-1}
                        i, i+1, diff_const*(-(1-C(i))-R(i)+v_2*(+3*(1-C(i))*C_diff-3*(1-C(i))^2)+v_3*(9*(1-C(i))^2*C_diff-5*(1-C(i))^3))-adv_const*C_adv_plus;                    %Coefficient terms of R_{i+1}
                        
                        i, i+row-1, -diff_const*R(i);                   %Coefficient terms of L_{i-1}
                        i, i+row, 2*diff_const*R(i)-react_const;        %Coefficient terms of L_{i}
                        i, i+row+1, -diff_const*R(i);                   %Coefficient terms of L_{i+1}
                        
                        i, i+2*row-1, -diff_const*R(i);                  %Coefficient terms of U_{i-1}
                        i, i+2*row, 2*diff_const*R(i)-react_const;       %Coefficient terms of U_{i,j}
                        i, i+2*row+1, -diff_const*R(i);                  %Coefficient terms of U_{i+1}
                        
                        i, i+3*row-1, -diff_const*R(i);                  %Coefficient terms of D_{i-1}
                        i, i+3*row, 2*diff_const*R(i)-react_const;       %Coefficient terms of D_{i,j}
                        i, i+3*row+1, -diff_const*R(i);                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            % Second quorter
            for i=row+1:2*row 
                ii=i-row;
                
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    
                    %Define the finite difference of the total density
                    %C_diff
                    C_diff=(C(ii+1)-C(ii-1))/2;
                    
                    if C(ii+1)~=0
                        C_adv_plus=(1-C(ii+1))*((1-C(ii+1))^v-1)/(C(ii+1));
                    else
                        C_adv_plus=0;
                    end
                    if C(ii-1)~=0
                        C_adv_minus=(1-C(ii-1))*((1-C(ii-1))^v-1)/(C(ii-1));
                    else
                        C_adv_minus=0;
                    end
                    
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_const*(2*(1-C(ii))+2*L(ii)+v_2*6*(1-C(ii))^2+v_3*(10*(1-C(ii))^3-6*(1-C(ii))*C_diff^2+3*(1-C(ii))^2*(C(ii+1)-2*C(ii)+C(ii-1))))+3*react_const;       %Coefficient terms of R_{i}
                        i, i-1, diff_const*(-(1-C(ii))-L(ii)+v_2*(-3*(1-C(ii))*C_diff-3*(1-C(ii))^2)+v_3*(-9*(1-C(ii))^2*C_diff-5*(1-C(ii))^3))-adv_const*C_adv_minus;                     %Coefficient terms of R_{i-1}
                        i, i+1, diff_const*(-(1-C(ii))-L(ii)+v_2*(+3*(1-C(ii))*C_diff-3*(1-C(ii))^2)+v_3*(9*(1-C(ii))^2*C_diff-5*(1-C(ii))^3))+adv_const*C_adv_plus;                       %Coefficient terms of R_{i+1}
                        
                        i, i-row-1, -diff_const*L(ii);                   %Coefficient terms of L_{i-1}
                        i, i-row, 2*diff_const*L(ii)-react_const;        %Coefficient terms of L_{i}
                        i, i-row+1, -diff_const*L(ii);                   %Coefficient terms of L_{i+1}
                        
                        i, i+1*row-1, -diff_const*L(ii);                  %Coefficient terms of U_{i-1}
                        i, i+1*row, 2*diff_const*L(ii)-react_const;       %Coefficient terms of U_{i,j}
                        i, i+1*row+1, -diff_const*L(ii);                  %Coefficient terms of U_{i+1}
                        
                        i, i+2*row-1, -diff_const*L(ii);                  %Coefficient terms of D_{i-1}
                        i, i+2*row, 2*diff_const*L(ii)-react_const;       %Coefficient terms of D_{i,j}
                        i, i+2*row+1, -diff_const*L(ii);                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            % Third quorter
            for i=2*row+1:3*row 
                ii=i-2*row;
                
               %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    
                    %Define the finite difference of the total density
                    %C_diff
                    C_diff=(C(ii+1)-C(ii-1))/2;
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_const*(2*(1-C(ii))+2*U(ii)+v_2*6*(1-C(ii))^2+v_3*(10*(1-C(ii))^3-6*(1-C(ii))*C_diff^2+3*(1-C(ii))^2*(C(ii+1)-2*C(ii)+C(ii-1))))+3*react_const;      %Coefficient terms of R_{i}
                        i, i-1, diff_const*(-(1-C(ii))-U(ii)+v_2*(-3*(1-C(ii))*C_diff-3*(1-C(ii))^2)+v_3*(-9*(1-C(ii))^2*C_diff-5*(1-C(ii))^3));                                          %Coefficient terms of R_{i-1}
                        i, i+1, diff_const*(-(1-C(ii))-U(ii)+v_2*(+3*(1-C(ii))*C_diff-3*(1-C(ii))^2)+v_3*(9*(1-C(ii))^2*C_diff-5*(1-C(ii))^3));                                           %Coefficient terms of R_{i+1}
                        
                        i, i-2*row-1, -diff_const*U(ii);                   %Coefficient terms of L_{i-1}
                        i, i-2*row, 2*diff_const*U(ii)-react_const;        %Coefficient terms of L_{i}
                        i, i-2*row+1, -diff_const*U(ii);                   %Coefficient terms of L_{i+1}
                        
                        i, i-1*row-1, -diff_const*U(ii);                  %Coefficient terms of U_{i-1}
                        i, i-1*row, 2*diff_const*U(ii)-react_const;       %Coefficient terms of U_{i,j}
                        i, i-1*row+1, -diff_const*U(ii);                  %Coefficient terms of U_{i+1}
                        
                        i, i+1*row-1, -diff_const*U(ii);                  %Coefficient terms of D_{i-1}
                        i, i+1*row, 2*diff_const*U(ii)-react_const;       %Coefficient terms of D_{i,j}
                        i, i+1*row+1, -diff_const*U(ii);                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                    %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            % Fourth quorter
            for i=3*row+1:4*row 
                ii=i-3*row;
                
                %Indeces corresponding to non-boundary sites
                if mod(i,row)>1
                    
                    %Define the finite difference of the total density
                    %C_diff
                    C_diff=(C(ii+1)-C(ii-1))/2;
                    coord(counter+1:counter+12,:)=...
                        [i, i, 1+diff_const*(2*(1-C(ii))+2*D(ii)+v_2*6*(1-C(ii))^2+v_3*(10*(1-C(ii))^3-6*(1-C(ii))*C_diff^2+3*(1-C(ii))^2*(C(ii+1)-2*C(ii)+C(ii-1))))+3*react_const;       %Coefficient terms of R_{i}
                        i, i-1, diff_const*(-(1-C(ii))-D(ii)+v_2*(-3*(1-C(ii))*C_diff-3*(1-C(ii))^2)+v_3*(-9*(1-C(ii))^2*C_diff-5*(1-C(ii))^3));                                           %Coefficient terms of R_{i-1}
                        i, i+1, diff_const*(-(1-C(ii))-D(ii)+v_2*(+3*(1-C(ii))*C_diff-3*(1-C(ii))^2)+v_3*(9*(1-C(ii))^2*C_diff-5*(1-C(ii))^3));                                            %Coefficient terms of R_{i+1}
                        
                        i, i-3*row-1, -diff_const*D(ii);                   %Coefficient terms of L_{i-1}
                        i, i-3*row, 2*diff_const*D(ii)-react_const;        %Coefficient terms of L_{i}
                        i, i-3*row+1, -diff_const*D(ii);                   %Coefficient terms of L_{i+1}
                        
                        i, i-2*row-1, -diff_const*D(ii);                  %Coefficient terms of U_{i-1}
                        i, i-2*row, 2*diff_const*D(ii)-react_const;       %Coefficient terms of U_{i,j}
                        i, i-2*row+1, -diff_const*D(ii);                  %Coefficient terms of U_{i+1}
                        
                        i, i-1*row-1, -diff_const*D(ii);                  %Coefficient terms of D_{i-1}
                        i, i-1*row, 2*diff_const*D(ii)-react_const;       %Coefficient terms of D_{i,j}
                        i, i-1*row+1, -diff_const*D(ii);                  %Coefficient terms of D_{i+1}
                        ];
                    counter=counter+12;
                    
                   %Indeces corresponding to boundary sites
                else
                    coord(counter+1,:)=[i, i, 1];
                    counter=counter+1;
                end % if
            end %for
            
            
            
            %% Solution of the linear system
            
            % Split the solution for the different polarisations
            R=u(1:row)';
            L=u(row+1:2*row)';
            U=u(2*row+1:3*row)';
            D=u(3*row+1:4*row)';
            
            %Define the vector b of the linear system
            %We force certain elements of the vector b to impose the boundary condition
            b=[0;R(2:end-1);0;0;L(2:end-1);0;0;U(2:end-1);0;0;D(2:end-1);0]; %                    
            
            %Build the sparse matrix A
            A=sparse(coord(:,1),coord(:,2),coord(:,3),size(u_IC,2),size(u_IC,2));
            
            %Solve the linear system
            u_plus=A\b;
            u_plus=u_plus';
            %Update auxiliary vector
            u_temp=u_plus;
            
            %Update number of iteration
            iter=iter+1;
            
        end %while
        %%
        
        %Display an error message if the tolerance has not been reached
        if iter==max_iter
            error('Picard Iteration reaches maximum iteration number');
        end % if
        
        % Update the Backward Euler time iteration
        u=u_plus;
        
    end % for (TIME LOOP)
       
end


% Split the solution for the different polarisations
R_n=u(1:row);
L_n=u(row+1:2*row);
U_n=u(2*row+1:3*row);
D_n=u(3*row+1:end);

%% Mass balance
fprintf('Mass balance: %2.2f percentege \n',(1-sum(u)/sum(u_IC))*100)

end

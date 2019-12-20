%
% Flow past NACA 0012
% Author: Shion Matsumoto
% Course: AEM 5253
%   Date: 12/13/2019

clear; clc; close all


%% MESHING
% Read mesh file
file = 'mesh_cooper_200x200.txt';
fid = fopen(file,'r');
formatSpec = '%d, %d, %f, %f \n';
szA = [4 Inf];
A = fscanf(fid,formatSpec,szA)';
fclose(fid);

% Parse data
xi = A(:,1);
eta = A(:,2);
x = A(:,3);
y = A(:,4);

% Rearrange data into matrix
% Rows - constant eta
% Columns - constant xi
xi = reshape(xi,200,[]);
eta = reshape(eta,200,[]);
x = reshape(x,200,[]);
y = reshape(y,200,[]);

% Plot mesh
figure
plot(x,y,'-k'); hold on;    % plot constant xi
plot(x',y','-k')            % plot constant eta
title('NACA0012 mesh in x-y plane','FontSize',18)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
axis equal tight

% Calculate Jacobian of grid transformation numerically
dxi = 1;
deta = 1;
n_xi = 200;
n_eta = 200;
metric_num = metric_calc_num(x,y,dxi,deta);

% Plot Jacobian (x-y plane)
plot_jacobian(x,y,xi,eta,metric_num.J,0);



%% SOLVE FOR FLOW PAST AIRFOIL
% Mach numbers for flow
mach = 0.85;

% Fluid properties
rho_inf = 1;
gamma = 1.4;
p_inf = 1 / gamma;
v_inf = 0;

% BC info
BC.rho_inf = rho_inf;
BC.gamma = gamma;
BC.p_inf = 1 / gamma;
BC.v_inf = v_inf;

% Physical domain ICs
rho0 = rho_inf * ones(n_eta,n_xi);
gamma = gamma * ones(n_eta,n_xi);
p0 = p_inf * ones(n_eta,n_xi);
v0 = zeros(n_eta,n_xi);

% Computational temporal and domain discretization
D.dt = 1e-4;    % time step
D.NT = 20000;   % number of time steps
D.dxi = dxi;
D.deta = deta;
D.N = n_xi;
D.M = n_eta;

% Convergence criteria (relative tolerance)
rtol = 1e-7;

% Initial free-stream velocity (x-dir)
u0 = mach * ones(n_eta,n_xi);

% Computational domain ICs
IC = calc_ic_metric(rho0,u0,v0,gamma,metric_num);

% BCs
BC.u_inf = mach;
BC.rhoE_inf = 1 / (BC.gamma * (BC.gamma - 1)) + (BC.u_inf ^ 2 + BC.v_inf ^ 2) / 2;

% Solve for UBAR

[UBAR,res] = solve_euler(IC,BC,D,metric_num,rtol,0);
semilogy(1:length(res),res); hold on

% Convert to Cartesian fields
U = UBAR .* metric_num.J;

rho = U(:,:,1);
u = U(:,:,2) ./ rho;
v = U(:,:,3) ./ rho;

figure
vel = sqrt(u .^ 2 + v .^ 2);
contourf(x,y,vel); hold on
%     quiver(x,y,u,v,0.5)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
title(['Steady-state velocity magnitude at M=' num2str(mach)],...
       'FontSize',18)
colorbar




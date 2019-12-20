% Name: Shion Matsumoto
% Date: 12/3/2019
% 
% AEM 5253 PROJECT

clear; clc; close all;


%% MESHING
% Computational domain discretization
n_xi = 100;
n_eta = 100;

% OPT = 0 (Whole domain)
% OPT = 1 (Top half domain)
OPT = 1;

switch OPT
    case 0  % Whole domain
        xi = linspace(0,2*pi,n_xi);
        eta = linspace(atanh(1/4),3,n_eta);
    case 1  % Half domain
        xi = linspace(0,pi,n_xi);
        eta = linspace(atanh(1/4),3,n_eta);
end

% Calculate step size
dxi = xi(2) - xi(1);
deta = eta(2) - eta(1);

% Cartesian coordinates of ellipse (4-to-1)
plot_mesh = 0;
plot_bc = 1;
[x,y] = create_mesh(xi,eta,plot_mesh,plot_bc);


%% METRIC TERMS (ANALYTICAL)
% Create meshgrid of eta and xi
[xixi,ee] = meshgrid(xi,eta);

% Calculate metric terms using analytic formula
metric_ana = metric_calc_ana(xixi,ee);

% Plot Jacobian
plot_jacobian(x,y,xi,eta,metric_ana.J,1);


%% METRIC TERMS (NUMERICAL)
% Calculate Jacobian of grid transformation numerically
metric_num = metric_calc_num(x,y,dxi,deta);

% Plot Jacobian (x-y plane)
plot_jacobian(x,y,xi,eta,metric_num.J,0);

% Compare analytical and numerical Jacobian
% J_err = sqrt((metric_ana.J - metric_num.J) .^ 2);
% 
% figure
% subplot(2,1,1)
% contourf(x,y,J_err); hold on
% title('Norm of numerical Jacobian error (x-y plane)','FontSize',14)
% xlabel('x','FontSize',14)
% ylabel('y','FontSize',14)
% colorbar
% 
% subplot(2,1,2)
% contourf(xi,eta,J_err); hold on
% title('Norm of numerical Jacobian error (\xi-\eta plane)','FontSize',14)
% xlabel('\xi','FontSize',14)
% ylabel('\eta','FontSize',14)
% colorbar


%% SOLVE FOR VARIOUS MACH NUMBERS
% Mach numbers for flow
mach = [0.4, 0.85, 1.5, 6];
mach = 0.4;

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
rtol = 1e-6;

% Solve U for each Mach number
figure
leg = {};
U = {};
for idx = 1:length(mach)
    
    % Initial free-stream velocity (x-dir)
    u0 = mach(idx) * ones(n_eta,n_xi);
    
    % Computational domain ICs
    IC = calc_ic_metric(rho0,u0,v0,gamma,metric_ana);
    
    % BCs
    BC.u_inf = mach(idx);
    BC.rhoE_inf = 1 / (BC.gamma * (BC.gamma - 1)) + (BC.u_inf ^ 2 + BC.v_inf ^ 2) / 2;
    
    % Solve for UBAR
    [UBAR,res] = solve_euler(IC,BC,D,metric_ana,rtol,OPT);
    semilogy(1:length(res),res); hold on
    leg{idx} = ['M=' num2str(mach(idx))];

    % Convert to Cartesian fields
    U{idx} = UBAR .* metric_ana.J;
end

% Residual plot
grid on
xlabel('Iteration No.')
ylabel('Residual')
legend(leg)
title('Residual Plot')

% Plot velocity fields
for idx = 1:length(mach)
    rho = U{idx}(:,:,1);
    u = U{idx}(:,:,2) ./ rho;
    v = U{idx}(:,:,3) ./ rho;

    figure
    vel = sqrt(u .^ 2 + v .^ 2);
    contourf(x,y,vel); hold on
%     quiver(x,y,u,v,0.5)
    xlabel('x','FontSize',14)
    ylabel('y','FontSize',14)
    title(['Steady-state velocity magnitude at M=' num2str(mach(idx))],...
           'FontSize',18)
    colorbar

    sonic = vel > 1;
    if any(any(sonic))
        figure
        contourf(x,y,sonic,[0 1])
        xlabel('x','FontSize',14)
        ylabel('y','FontSize',14)
        title(['Subsonic and supersonic regions for M=' num2str(mach(idx))],...
               'FontSize',18)
    end
end


% %% DRAG COEFFICIENT CALCULATION
% % Calculate drag force and drag coefficient for each Mach number
% 
% % Integrate x-direction of pressure over surface of ellipse to find drag
% % force
% 
% C_D = [];
% for idx = 1:length(mach)
%     % Get density, pressure, and energy
%     rho = U{idx}(:,:,1);
%     u = U{idx}(:,:,2) ./ rho;
%     v = U{idx}(:,:,3) ./ rho;
%     rhoE = U{idx}(:,:,4);
%     
%     % Calculate pressure field
%     p = (gamma - 1) .* (rhoE - (rho / 2) .* (u .^ 2 + v .^ 2));
%     
%     % Calculate drag coefficient
%     C_D(idx) = calc_drag(p(1,:),x(1,:),y(1,:),rho_inf,mach(idx));
% end
% 
% % Plot C_D as a function of Mach 0.4 to 1.5
% figure
% plot(mach,C_D,'-*'); hold on; grid on
% xlabel('Mach Number (M)','FontSize',14)
% ylabel('Drag Coefficient (C_D)','FontSize',14)
% title('Mach number vs. Drag coefficient','FontSize',18)


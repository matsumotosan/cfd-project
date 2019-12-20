function [rho,rhoE,u,v,UBAR] = enter_bc(UBAR,BC,metric,OPT)
%CALC_BC_METRIC Calculate boundary conditions for FBAR and GBAR
%   Boundary conditions:
%   1) Freestream edge (xi=0,2pi, eta=3): free stream velocity
%
%   2) Ellipse surface (xi=0 to 2pi, eta=0.25): free slip
%   - (du_tilda/deta = 0)
%   - v_tilda = 0
%

% Parse BC inputs                                 
gamma = BC.gamma;
rho_inf = BC.rho_inf;
u_inf = BC.u_inf;
v_inf = BC.v_inf;
rhoE_inf = BC.rhoE_inf;

% Calculate rho,u,v,E in x- and y- coordinates from UBAR
rho = UBAR(:,:,1) .* metric.J;
u = (UBAR(:,:,2) .* metric.J) ./ rho;
v = (UBAR(:,:,3) .* metric.J) ./ rho;
rhoE = UBAR(:,:,4) .* metric.J;

% Calculate contravariant velocities
u_tilda = metric.xi_x .* u + metric.xi_y .* v;
v_tilda = metric.eta_x .* u + metric.eta_y .* v;

switch OPT
    case 0  % Whole domain
        % Free-slip & no-flux boundary at ellipse surface
        u_tilda(1,:) = u_tilda(2,:);    % free-slip (du_tilda/deta=0)
        v_tilda(1,:) = 0;               % no flux (v_tilda=0)
        rho(1,:) = rho(2,:);            % (drho/deta=0)
        rhoE(1,:) = rhoE(2,:);          % (dE/deta=0)

        % Periodic boundary condition at xi=0 and xi=2pi
        u_tilda(:,end) = u_tilda(:,1);
        v_tilda(:,end) = v_tilda(:,1);
        rho(:,end) = rho(:,1);
        rhoE(:,end) = rhoE(:,1);

        % Convert back to physical velocities
        % My formula
        % u = (u_tilda ./ metric.xi_y - v_tilda ./ metric.eta_y) ./ ...
        %     (-metric.eta_x ./ metric.eta_y + metric.xi_x ./ metric.xi_y);
        % v = (u_tilda ./ metric.xi_x - v_tilda ./ metric.eta_x) ./ ...
        %     (-metric.eta_y ./ metric.eta_x + metric.xi_y ./ metric.xi_x);

        % Other formula to convert back to physical velocities
        u = (metric.eta_y .* u_tilda - metric.xi_y .* v_tilda) ./ metric.J;
        v = (metric.xi_x .* v_tilda - metric.eta_x .* u_tilda) ./ metric.J;

        % Far-field free values
        rho(end,:) = rho_inf;
        u(end,:) = u_inf;
        v(end,:) = v_inf;
        rhoE(end,:) = rhoE_inf;
        
    case 1  % Half domain
        % Free-slip & no-flux boundary at ellipse surface
        u_tilda(1,:) = u_tilda(2,:);    % free-slip (du_tilda/deta=0)
        v_tilda(1,:) = 0;               % no flux (v_tilda=0)
        rho(1,:) = rho(2,:);            % Neumann BC (drho/deta=0)
        rhoE(1,:) = rhoE(2,:);          % Neumann BC (drhoE/deta=0)
        
        % Right side of ellipse (symmetric/Neumann BC)
        u_tilda(:,1) = 0;
        v_tilda(:,1) = v_tilda(:,2);
        rho(:,1) = rho(:,2);
        rhoE(:,1) = rhoE(:,2);
        
        % Left side of ellipse (symmetric/Neumann BC)
        u_tilda(:,end) = 0;
        v_tilda(:,end) = v_tilda(:,end-1);
        rho(:,end) = rho(:,end-1);
        rhoE(:,end) = rhoE(:,end-1);
        
%         u = (u_tilda ./ metric.xi_y - v_tilda ./ metric.eta_y) ./ ...
%             (-metric.eta_x ./ metric.eta_y + metric.xi_x ./ metric.xi_y);
%         v = (u_tilda ./ metric.xi_x - v_tilda ./ metric.eta_x) ./ ...
%             (-metric.eta_y ./ metric.eta_x + metric.xi_y ./ metric.xi_x);

        % Convert back to physical velocities (Alice's formula)
        u = (metric.eta_y .* u_tilda - metric.xi_y .* v_tilda) ./ metric.J;
        v = (metric.xi_x .* v_tilda - metric.eta_x .* u_tilda) ./ metric.J;

        % Far-field values (Dirichlet BC)
        rho(end,:) = rho_inf;
        u(end,:) = u_inf;
        v(end,:) = v_inf;
        rhoE(end,:) = rhoE_inf;
end

% Update UBAR
UBAR = cat(3,rho,rho .* u,rho .* v,rhoE) ./ metric.J;

end


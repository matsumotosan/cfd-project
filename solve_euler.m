function [UBAR,res] = solve_euler(IC,BC,D,metric,rtol,OPT)
%SOLVE_EULER Solve Euler equations in generalized curvilinear
% coordinates using a modified version of a scheme from Steger 
% and Warming (Journal of Computational Physics, 1981). 
% Upwind differencing in space and Euler stepping (first order)
% in time. 

% Initialize solution matrix and enter initial conditions
UBAR = IC.UBAR;
gamma = IC.gamma;
[N,M,~] = size(UBAR);
res = zeros(D.NT,1);

% Solve for U numerically
for I = 1:D.NT
    if ~isreal(UBAR)
        disp('UBAR is complex.')
        break
    end
    
    % Enter BCs into UBAR and return velocities in physical domain
    [rho,rhoE,u,v,UBAR] = enter_bc(UBAR,BC,metric,OPT);
    
    % Flux vector splitting (Steger & Warming, 1981)
    [Fp,Fn,Gp,Gn] = flux_split(gamma,rho,rhoE,u,v,metric);
    
    % Upwind difference to calculate fluxes
    [dF_dxi,dG_deta] = upwind_diff(Fp,Fn,Gp,Gn,D.dxi,D.deta);
    
    % Solve for next time step
    UBAR_new = UBAR - D.dt * (dF_dxi + dG_deta);
    
    % Calculate residual
    res(I) = calc_res(UBAR_new,UBAR,N,M);
    
    % Check if converged
    if res(I) < rtol
        idx = find(res == 0,1,'first');
        res = res(1:idx);
        UBAR = UBAR_new;
        fprintf('SOLUTION CONVERGED.\n')
        return
    end
    
    % Update UBAR
    UBAR = UBAR_new;
    
%     if mod(I,100) == 1
%         U = UBAR .* metric.J;
%     end
    
    % Check for NaN in solution
    if any(isnan(UBAR(:)))
        fprintf('SOLUTION HAS NAN. ITERATION:%d\n',I);
        return
    end
end

% Calculate residual if solution did not converge
fprintf('SOLUTION DID NOT CONVERGE.\nResidual:%.4f\n',res(end));

end


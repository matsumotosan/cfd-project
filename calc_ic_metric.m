function IC = calc_ic_metric(rho0,u0,v0,gamma,metric)
%CALC_IC_METRIC Calculate initial conditions for UBAR.
%   Returns initial conditions in a struct.

p0 = 1 ./ gamma;
rhoE0 = 1 ./ (gamma .* (gamma - 1)) + (u0 .^ 2) ./ 2;
UBAR = cat(3,rho0,rho0 .* u0,rho0 .* v0,rhoE0) ./ metric.J;
c = sqrt(gamma .* p0 ./ rho0);

% Store initial conditions in struct
IC.UBAR = UBAR;
IC.gamma = gamma;
IC.c = c;
IC.p0 = p0;

end


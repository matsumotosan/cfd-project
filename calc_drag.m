function C_D = calc_drag(p,x,y,rho_inf,U_inf)
%CALC_PRES Calculate drag coefficient
%   Numerically integrate x-component of pressure around ellipse and
%   calculate drag force F_D.

% Calculate surface normal unit vector
n = calc_normal(x,y);

% Numerically integrate x-component of pressure over ellipse surface
dl = calc_len(x,y);

% Find average pressure between consecutive grid points
p = 0.5 * (p(1:end-1) + p(2:end));

% Numerically integrate pressure in x-direction around ellipse for F_D
e1 = repmat([1,0],[length(n),1]);
F_D = sum(-p(:) .* dot(e1,n,2) .* dl(:));

% Calculate drag coefficient C_D
A = sum(dl);
C_D = 2 * F_D / (rho_inf * U_inf ^ 2 * A);

end


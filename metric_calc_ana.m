function metric_ana = metric_calc_ana(xi,eta)
%METRIC_CALC Calculate metric terms for x and y in terms of 
%   eta and xi analytically.

% Calculate metric terms for x and y
x_xi = -cosh(eta) .* sin(xi);
x_eta = sinh(eta) .* cos(xi);

y_xi = sinh(eta) .* cos(xi);
y_eta = cosh(eta) .* sin(xi);

% Calculate Jacobian of grid transformation
J = 1 ./ (x_xi .* y_eta - y_xi .* x_eta);

% Return metric terms in struct
metric_ana.x_xi = x_xi;
metric_ana.x_eta = x_eta;
metric_ana.y_xi = y_xi;
metric_ana.y_eta = y_eta;
metric_ana.xi_x = J .* y_eta;
metric_ana.xi_y = -J .* x_eta;
metric_ana.eta_x = -J .* y_xi;
metric_ana.eta_y = J .* x_xi;
metric_ana.J = J;

end


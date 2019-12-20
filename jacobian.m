function metric = jacobian(x,y,xi,eta)
%JACOBIAN Numerically calculate Jacobian of grid transformation.

% Calculate metric terms numerically
[x_xi,x_eta,y_xi,y_eta] = metric_calc_num(x,y,xi,eta);

% Jacobian of grid transformation
J = 1 ./ (x_xi .* y_eta - y_xi .* x_eta);

% Return metric terms in struct
metric.x_xi = x_xi;
metric.x_eta = x_eta;
metric.y_xi = y_xi;
metric.y_eta = y_eta;
metric.xi_x = J .* y_eta;
metric.xi_y = -J .* x_eta;
metric.eta_x = -J .* y_xi;
metric.eta_y = J .* x_xi;
metric.J = J;

end


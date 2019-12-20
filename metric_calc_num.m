function metric_num = metric_calc_num(x,y,dxi,deta)
%METRIC_CALC_NUM Numerical calculation of metric terms for
%   interior and boundary grid points. Calculate following 
%   metric terms numerically:
%       - x_xi
%       - x_eta
%       - y_xi
%       - y_eta
%
%   Function assumes that d(eta) and d(xi) are constant. 

% Initialize metric terms
x_xi = zeros(size(x));
x_eta = x_xi;
y_xi = x_xi;
y_eta = x_xi;

% Metric terms for interior nodes
x_xi(:,2:end-1) = (x(:,3:end) - x(:,1:end-2)) / (2 * dxi);
x_eta(2:end-1,:) = (x(3:end,:) - x(1:end-2,:)) / (2 * deta);

y_xi(:,2:end-1) = (y(:,3:end) - y(:,1:end-2)) / (2 * dxi);
y_eta(2:end-1,:) = (y(3:end,:) - y(1:end-2,:)) / (2 * deta);

% Metric terms for left boundary nodes
x_xi(:,1) = (-3 * x(:,1) + 4 * x(:,2) - x(:,3)) / (2 * dxi);
y_xi(:,1) = (-3 * y(:,1) + 4 * y(:,2) - y(:,3)) / (2 * dxi);

% Metric terms for right boundary nodes
x_xi(:,end) = (3 * x(:,end) - 4 * x(:,end-1) + x(:,end-2)) / (2 * dxi);
y_xi(:,end) = (3 * y(:,end) - 4 * y(:,end-1) + y(:,end-2)) / (2 * dxi);

% Metric terms for bottom boundary nodes
x_eta(1,:) = (-3 * x(1,:) + 4 * x(2,:) - x(3,:)) / (2 * deta);
y_eta(1,:) = (-3 * y(1,:) + 4 * y(2,:) - y(3,:)) / (2 * deta);

% Metric terms for top boundary nodes
x_eta(end,:) = (3 * x(end,:) - 4 * x(end-1,:) + x(end-2,:)) / (2 * deta);
y_eta(end,:) = (3 * y(end,:) - 4 * y(end-1,:) + y(end-2,:)) / (2 * deta);

% Jacobian of grid transformation
J = 1 ./ (x_xi .* y_eta - y_xi .* x_eta);

% Return metric terms in struct
metric_num.x_xi = x_xi;
metric_num.x_eta = x_eta;
metric_num.y_xi = y_xi;
metric_num.y_eta = y_eta;
metric_num.xi_x = J .* y_eta;
metric_num.xi_y = -J .* x_eta;
metric_num.eta_x = -J .* y_xi;
metric_num.eta_y = J .* x_xi;
metric_num.J = J;

end
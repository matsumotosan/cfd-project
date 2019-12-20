function plot_velocity(u,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_xi = 100;
n_eta = 50;

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
[x,y] = create_mesh(xi,eta,1);


figure
surf(x,y,sqrt(u .^ 2 + v .^ 2)); hold on
% quiver(x,y,u,v)
xlabel('x')
ylabel('y')
title('Velocity contour plot with vectors')
colorbar 

end


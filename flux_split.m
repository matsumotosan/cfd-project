 function [Fp,Fn,Gp,Gn] = flux_split(gamma,rho,rhoE,u,v,metric)
%FLUX_SPLIT Calculate flux vectors FP,FN,GP,GN for flux splitting

% Calculate speed of sound
p = (gamma - 1) .* (rhoE - (rho / 2) .* (u .^ 2 + v .^ 2));
c = sqrt(gamma .* p ./ rho);

if any(isnan(c(:)))
    fprintf('C has NaN.\n');
elseif any(~isreal(c(:)))
    fprintf('C is complex.\n');
end

% Find eigenvalues
F_lambda = find_eigen(metric.xi_x,metric.xi_y,u,v,c);
G_lambda = find_eigen(metric.eta_x,metric.eta_y,u,v,c);

% Split eigenvalues
a = 1e-7;
[F_lambda_plus,F_lambda_neg] = eigen_split(F_lambda,a);
[G_lambda_plus,G_lambda_neg] = eigen_split(G_lambda,a);

% Flux-splitter
Fp = flux_gen(F_lambda_plus,rho,gamma,c,metric.xi_x,metric.xi_y,u,v) ./ metric.J.^ 2;
Fn = flux_gen(F_lambda_neg,rho,gamma,c,metric.xi_x,metric.xi_y,u,v) ./ metric.J .^ 2;
Gp = flux_gen(G_lambda_plus,rho,gamma,c,metric.eta_x,metric.eta_y,u,v) ./ metric.J .^ 2;
Gn = flux_gen(G_lambda_neg,rho,gamma,c,metric.eta_x,metric.eta_y,u,v) ./ metric.J .^ 2;

% Fp = (metric.xi_x .* Fp + metric.xi_y .* Gp) ./ metric.J;
% Fn = (metric.xi_x .* Fn + metric.xi_y .* Gn) ./ metric.J;
% Gp = (metric.eta_x .* Fp + metric.eta_y .* Gp) ./ metric.J;
% Gn = (metric.eta_x .* Fn + metric.eta_y .* Gn) ./ metric.J;

end


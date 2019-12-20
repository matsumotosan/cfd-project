function F = flux_gen(lambda,rho,gamma,c,k1,k2,u,v)
%FLUX_GENERAL Calculate generalized flux vector F for scheme
%   by Steger and Warming (Journal of Computational Physics, 1981).
%   

% Calculate k1 and k2
k1_til = k1 ./ sqrt(k1 .^ 2 + k2 .^ 2);
k2_til = k2 ./ sqrt(k1 .^ 2 + k2 .^ 2);

l1 = lambda(:,:,1);
l3 = lambda(:,:,2);
l4 = lambda(:,:,3);

% Calculate flux vector components
F1 = 2 * (gamma - 1) .* l1 + l3 + l4;
F2 = 2 * (gamma - 1) .* l1 .* u + l3 .* (u + c .* k1_til) + ...
     l4 .* (u - c .* k1_til);
F3 = 2 * (gamma - 1) .* l1 .* v + l3 .* (v + c .* k2_til) + ...
     l4 .* (v - c .* k2_til);
W = ((3 - gamma) .* (l3 + l4) .* c .^ 2) ./ (2 * (gamma - 1));
F4 = (gamma - 1) .* l1 .* (u .^ 2 + v .^ 2) + (l3 / 2) .* ...
     ((u + c .* k1_til) .^ 2 + (v + c .* k2_til) .^ 2) + ...
     (l4 / 2) .* ((u - c .* k1_til) .^ 2 + (v - c .* k2_til) .^ 2) + W;

% Return flux vector
F = (rho ./ (2 * gamma)) .* cat(3,F1,F2,F3,F4);

end


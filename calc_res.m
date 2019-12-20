function res = calc_res(U,U_new,N,M)
%CALC_RES Summary of this function goes here
%   Detailed explanation goes here

% L2 norm for convergence criteria
res = sqrt(sum((U - U_new) .^ 2,'all')) / (N * M);

end


function [dF_dxi,dG_deta] = upwind_diff(Fp,Fn,Gp,Gn,dxi,deta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initialize matrices
dFp = zeros(size(Fp));
dFn = dFp;
dGp = dFp;
dGn = dFp;

% Upwind differencing
% dFp(:,2:end-1,:) = Fp(:,2:end-1,:) - Fp(:,1:end-2,:);
% dFn(:,2:end-1,:) = Fn(:,3:end,:) - Fn(:,2:end-1,:);
% dGp(2:end-1,:,:) = Gp(2:end-1,:,:) - Gp(1:end-2,:,:);
% dGn(2:end-1,:,:) = Gn(3:end,:,:) - Gn(2:end-1,:,:);

dFp(2:end-1,2:end-1,:) = Fp(2:end-1,2:end-1,:) - Fp(2:end-1,1:end-2,:);
dFn(2:end-1,2:end-1,:) = Fn(2:end-1,3:end,:) - Fn(2:end-1,2:end-1,:);
dGp(2:end-1,2:end-1,:) = Gp(2:end-1,2:end-1,:) - Gp(1:end-2,2:end-1,:);
dGn(2:end-1,2:end-1,:) = Gn(3:end,2:end-1,:) - Gn(2:end-1,2:end-1,:);

% Calculate upwind differences for each flux term
dF_dxi = (dFp + dFn) / dxi;
dG_deta = (dGp + dGn) / deta;

end


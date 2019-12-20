function [lambda_plus,lambda_neg] = eigen_split(lambda,a)
%EIGENSPLIT Split eigenvalues into positive and negative without
%   corners at zero.
%

lambda_plus = (lambda + sqrt(lambda .^ 2 + a ^ 2)) / 2;
lambda_neg = (lambda - sqrt(lambda .^ 2 + a ^ 2)) / 2;

end


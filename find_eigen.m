function lambda = find_eigen(k1,k2,u,v,c)
%FINDEIGEN Summary of this function goes here
%   Detailed explanation goes here

l1 = k1 .* u + k2 .* v;
l3 = l1 + c .* sqrt(k1 .^ 2 + k2 .^ 2);
l4 = l1 - c .* sqrt(k1 .^ 2 + k2 .^ 2);

lambda = cat(3,l1,l3,l4);

end


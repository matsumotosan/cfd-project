function n = calc_normal(x,y)
%CALC_NORMAL Calculate surface normal unit vector and length of line
%elements

% Calculate vector pointing outward
dx = diff(x);
dy = diff(y);
len = sqrt(dx .^ 2 + dy .^ 2);
n = [dy(:)./ len(:),-dx(:) ./ len(:)];

end


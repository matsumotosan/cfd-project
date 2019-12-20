function dl = calc_len(x,y)
%CALC_LENGTH Calculate length of line elements between grid points on
%ellipse

dl = sqrt(diff(x) .^ 2 + diff(y) .^ 2);

end


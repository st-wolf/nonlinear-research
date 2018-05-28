function d = dissipation(xgrid)

dissip_stength = 500;

l = 0.2 * (xgrid(end) - xgrid(1));
L = xgrid(end);

d = zeros(1, length(xgrid));

sel = (xgrid > (L - l));
d(sel) = dissip_stength * ((xgrid(sel) - (L - l)) .^ 2);

sel = (xgrid < (-L + l));
d(sel) = dissip_stength * ((xgrid(sel) + (L - l)) .^ 2);

end


function  Y = resample(X0, Y0, X)
% Interpolate the solution (Y0) from one grid (X0) to another (X)

Y = zeros(size(X));

sel = find((X <= X0(end)) & (X >= X0(1)));
Y(sel) = spline(X0, Y0, X(sel));

sel = find((X > X0(end)) | (X < X0(1)));
if ~isempty(sel)
	Y(sel) = 0; %zero tails
end;

end

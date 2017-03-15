function x_root = newton( f, df, x0 )
%
% INPUT:
%	

eps = 1e-6;

while true
	x1 = x0 - f(x0) / df(x0);
	
	if abs(x1 - x0) < eps
		break
	else
		x0 = x1;
	end
end

x_root = x1;

end


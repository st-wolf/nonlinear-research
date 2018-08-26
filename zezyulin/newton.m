function x_root = newton( f, x0 )
%
% INPUT:
%	

MAX_ITERATION = 50;

eps = 1e-6;

df = @(x) (f(x + eps) - f(x - eps)) / (2 * eps);

iter = 0;
while true
	iter = iter + 1;
	
	x1 = x0 - f(x0) / df(x0);
	
	if abs(x1 - x0) < eps
		break
	else
		x0 = x1;
	end
	
	if iter > MAX_ITERATION
		fprintf('--> MAX_ITERATION\n')
		break
	end
end

x_root = x1;

end


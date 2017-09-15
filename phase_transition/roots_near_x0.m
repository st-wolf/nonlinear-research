function [ xleft, xright ] = roots_near_x0( f, x0 )
%
% INPUT:
%

% Step
step = abs(f(x0));

% Convergence
delta = 1e-12;

x_prev = x0;
while true
	x_next = x_prev + step;
	
	if sign(f(x_prev)) == sign(f(x_next))
		x_prev = x_next;
	else
		xright = dichotomy(f, x_prev, x_next, delta);
		break
	end
end

x_prev = x0;
while true
	x_next = x_prev - step;
	
	if sign(f(x_prev)) == sign(f(x_next))
		x_prev = x_next;
	else
		xleft = dichotomy(f, x_next, x_prev, delta);
		break
	end
end

end
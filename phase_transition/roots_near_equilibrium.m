function [ xleft, xright ] = roots_near_equilibrium( f )
%
% INPUT:
%

% Step
step = 0.01;

% Convergence
delta = 1e-12;

x0 = 0;
x1 = x0 + step;

while true
	if sign(f(x0)) == sign(f(x1))
		x0 = x1;
		x1 = x1 + step;
	else
		x_root = dichotomy(f, x0, x1, delta);
		break
	end
end

xright = abs(x_root);
xleft = -xright;

end
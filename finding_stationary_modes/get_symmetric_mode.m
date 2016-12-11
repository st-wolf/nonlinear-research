function [ X, Phi ] = get_symmetric_mode( params, C, xspan )
% Calculate a symmetric mode using the parameter C of the asymptotic 'asympt'
%
% INPUT:
%

xstart = xspan(1);

% From the asymtote of the solutions
asympt_params = @(x) asympt(params, C, x);

% Derivative
eps = 1e-12;

init = [ asympt_params(xstart)
		(asympt_params(xstart + eps) - asympt_params(xstart - eps)) / eps];

ode_params = @(x, f) ode(x, f, params);
[X, Phi] = RK4(ode_params, xspan, init, 1024);

X = [X -X(end:-1:1)];
Phi = [
	Phi(:, 1), Phi(:, 2);
	Phi(end:-1:1, 1), -Phi(end:-1:1, 2)
];

end
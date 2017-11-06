function [ X_approx, U_approx ] = get_another_alfimov_approximation( params, n, xspan )
% Equation: u_{xx} + (\mu - x^2) u + \sigma_1 \cos(\Omega x) u^3 = 0 
%
% INPUT:
%	params - [\mu \Omega sigma_1]
%	n - branch number: 0, 1, 2, 3 ...
%	x - grid, x = xstart:xstep:xend

% Unpacking
xstart = xspan(1); mu = params(1); Omega = params(2); sigma_1 = params(3);
mu_analog = 2*n + 1; % from \mu_n to \mu.

[X_approx, U] = get_mode_with_linear_counterpart('f_approx_solve', mu_analog, params, xstart);

% Approximation formula
% u(x) = \sqrt{\Omega} U(x) + (\sigma_1 / \sqrt{\Omega}) U^3(x) \cos(\Omega x)

% U_approx = ...
% 	sqrt(Omega) * U(:, 1);

U_approx = ...
	sqrt(Omega) * U(:, 1) + ...
	(sigma_1 / sqrt(Omega)) * (U(:, 1) .^ 3) .* cos(Omega * X_approx);

end


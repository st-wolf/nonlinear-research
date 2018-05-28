function u_approx = get_alfimov_approximation( params, n, x )
% Equation: u_{xx} + (\mu - x^2) u + \sigma_1 \cos(\Omega x) u^3 = 0 
%
% INPUT:
%	params - [\mu \Omega sigma_1]
%	n - branch number: 0, 1, 2, 3 ...
%	x - grid, x = xstart:xstep:xend

% Unpacking
xstep = x(2) - x(1); mu = params(1); Omega = params(2); sigma_1 = params(3);


mu_n = 2*n + 1; % from \mu_n to \mu. \Delta \mu = \mu_n - \mu
delta_mu = mu_n - mu;

Hn_solution = (1 / sqrt(2^n * factorial(n) * sqrt(pi))) * hermite(n, x) .* exp(-0.5 * (x .^ 2));
Hn_2nd_integral = simpson(Hn_solution .^ 2, xstep);
Hn_6th_integral = simpson(Hn_solution .^ 6, xstep);

U0 = ((2 * Hn_2nd_integral * delta_mu) / (3 * Hn_6th_integral)) ^ (1/4);

u_approx = ...
	U0 * sqrt(Omega) * Hn_solution + ...
	(U0 ^ 3) * (sigma_1 / sqrt(Omega)) * (Hn_solution .^ 3) .* cos(Omega * x);

% u_approx = ...
% 	U0 * sqrt(Omega) * Hn_solution;

end


function [ X, U ] = get_mode_with_linear_counterpart( mu_analog, params, xstart )
% For the equation: u_{xx} + (\mu - x^2) u + \sigma_1 \cos(\Omega x) u^3 = 0 
%
% INPUT:
%

mu_target = params(1); Omega = params(2); sigma_1 = params(3);
mu_aug = 0.05; % augmentation to bifurcate from zero solution
mu_start = mu_analog - mu_aug;

% Numer of the mode: 0, 1, 2, 3, 4, 5 ...
n = round((mu_analog - 1) / 2); % \mu = 2n + 1

if mod(n, 2) == 0
	% Even solution
	get_end = @(mex_solver_name, params, C, xspan) get_ux_end(mex_solver_name, params, C, xspan);
else
	% Odd solution
	get_end = @(mex_solver_name, params, C, xspan) get_u_end(mex_solver_name, params, C, xspan);
end

% Find solution in \mu_start and continue it on the parameter \mu
params = [mu_start Omega sigma_1];
get_end_params = @(c) get_end('f_sigma_solve', params, c, [xstart 0]);

% Find cspan for dichotomy solver
cstep = 0.1; cstart = cstep; cend = cstart + cstep;
end_cstart = get_end_params(cstart);
end_cend = get_end_params(cend);

while sign(end_cstart) == sign(end_cend)
	cend = cend + cstep;
	end_cend = get_end_params(cend);
end

% Dichotomy precision
eps = 1e-9;

cmode = dichotomy(get_end_params, cstart, cend, eps);

if mod(n, 2) == 0
	% Even solution
	[X, U] = get_symmetric_mode('f_sigma_solve', params, cmode, [xstrat 0]);
else
	% Odd solution
	[X, U] = get_antisymmetric_mode('f_sigma_solve', params, cmode, [xstart 0]);
end

plot(X, U);

% Continuation on the parameter \mu
mu_step = mu_aug; % why not!?

mu_intermediate_values = mu_start:(-mu_step):mu_target;
cmode_intermediate_values = zeros(1, length(mu_intermediate_values));
cmode_intermediate_values(1) = cmode;

for i = 2:length(mu_intermediate_values)
	fprintf('%i of %i\n', i, length(mu_intermediate_values))
	params = [mu_intermediate_values(i) Omega sigma_1];
	
	get_end_params = @(c) get_u_end('f_sigma_solve', params, c, [xstart 0]);
	cmode_intermediate_values(i) = newton(get_end_params, cmode_intermediate_values(i - 1));
	
	% Debug usage
	if mod(n, 2) == 0
		% Even solution
		[X, U] = get_symmetric_mode('f_sigma_solve', params, cmode_intermediate_values(i), [xstrat 0]);
	else
		% Odd solution
		[X, U] = get_antisymmetric_mode('f_sigma_solve', params, cmode_intermediate_values(i), [xstart 0]);
	end
	
	plot(X, U);
	pause()
end

end


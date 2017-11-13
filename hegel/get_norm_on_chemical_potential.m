function [ mu_intermediate_values, mode_norm, stability ] = ...
	get_norm_on_chemical_potential( mex_solver_name, mu_analog, params, xstart )
% Input arguments are the same as for @get_mode_with_linear_counterpart
% function. Nonlinear potential is P(x) = 1 + P_1 \cos(\Omega x) used in
% @get_spectrum function.
%
% Equation is u_{xx} + (\mu - x^2) u + (1 + P_1 \cos(\Omega x)) u^3 = 0
%
% INPUT:
%	mu_analog - linear counterpart value of chemical potential
%	params    - [mu_target Omega P1], where mu_target - final value of chemical
%		potential for the nonlinear mode
%
% OUTPUT:
%	mu_intermediate_values - chemical potential
%	mode_norm - norm of the solution
%	stability - vactor of 1s (stable) and 0s (unstable)
%

% -----------------------------------------------------------------------------
% Nonlinear potential
nonlinear_potential = @(params, x) 1 + params(3) * cos(params(2) * x);
% -----------------------------------------------------------------------------

mu_target = params(1); Omega = params(2); P1 = params(3);
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
params = [mu_start Omega P1];
get_end_params = @(c) get_end(mex_solver_name, params, c, [xstart 0]);

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
	[X, U] = get_symmetric_mode(mex_solver_name, params, cmode, [xstart 0]);
else
	% Odd solution
	[X, U] = get_antisymmetric_mode(mex_solver_name, params, cmode, [xstart 0]);
end

% Continuation on the parameter \mu
mu_step = mu_aug; % why not!?

mu_intermediate_values = mu_start:(-mu_step):mu_target;

if mu_intermediate_values(end) ~= mu_target
	mu_intermediate_values = [mu_intermediate_values mu_target];
end

cmode_intermediate_values = zeros(1, length(mu_intermediate_values));
cmode_intermediate_values(1) = cmode;

mode_norm = zeros(1, length(mu_intermediate_values));
stability = zeros(1, length(mu_intermediate_values));

% For the first mode
mode_norm(1) = get_norm(X, U);

% TODO: what are the parameters?
% The same as in @get_spectrum function
n_harmonics = 500;
stability(1) = is_stable(params, nonlinear_potential, X, U, n_harmonics);

for i = 2:length(mu_intermediate_values)
	fprintf('\tFrom linear analog: %i of %i\n', i, length(mu_intermediate_values))
	params = [mu_intermediate_values(i) Omega P1];
	
	get_end_params = @(c) get_end(mex_solver_name, params, c, [xstart 0]);
	cmode_intermediate_values(i) = newton(get_end_params, cmode_intermediate_values(i - 1));
	
	% Debug usage
	if mod(n, 2) == 0
		% Even solution
		[X, U] = get_symmetric_mode(mex_solver_name, params, cmode_intermediate_values(i), [xstart 0]);
	else
		% Odd solution
		[X, U] = get_antisymmetric_mode(mex_solver_name, params, cmode_intermediate_values(i), [xstart 0]);
	end
	
	% Computing norm and stability of the solution
	mode_norm(i) = get_norm(X, U);
	stability(i) = is_stable(params, nonlinear_potential, X, U, n_harmonics);
	
	% For debug purposes
	% plot(X, U);
	% pause()
end

end


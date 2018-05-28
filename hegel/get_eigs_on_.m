function [ all_mu, all_eigs ] = ...
	get_eigs_on_( mex_solver_name, mu_analog, params, xstart )

% -----------------------------------------------------------------------------
% Nonlinear potential
% nonlinear_potential = @(params, x) 1 + params(3) * cos(params(2) * x);
nonlinear_potential = @(params, x) params(3) * cos(params(2) * x);
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Imaginary window: [-12; +12]
get_from_window = @(x) x(abs(imag(x)) < 13);
% -----------------------------------------------------------------------------

mu_target = params(1); Omega = params(2); P1 = params(3);
mu_aug = 0.025; % augmentation to bifurcate from zero solution
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
% mu_step = mu_aug; % why not!?
mu_step = 0.025; % why not!?

mu_intermediate_values = mu_start:(-mu_step):mu_target;

if mu_intermediate_values(end) ~= mu_target
	mu_intermediate_values = [mu_intermediate_values mu_target];
end

cmode_intermediate_values = zeros(1, length(mu_intermediate_values));
cmode_intermediate_values(1) = cmode;

all_eigs = [];
all_mu = [];

% TODO: what are the parameters?
% The same as in @get_spectrum function
n_harmonics = 500;

eigenvalues = get_spectrum(params, nonlinear_potential, X, U, n_harmonics);
to_add = get_from_window(eigenvalues);

all_eigs = [all_eigs; to_add];
all_mu = [all_mu; mu_start * ones(length(to_add), 1)];

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
	
	eigenvalues = get_spectrum(params, nonlinear_potential, X, U, n_harmonics);
	to_add = get_from_window(eigenvalues);

	all_eigs = [all_eigs; to_add];
	all_mu = [all_mu; mu_intermediate_values(i) * ones(length(to_add), 1)];

	% plot_spectrum(params, to_add);
	% pause()
	
	% For debug purposes
	% plot(X, U);
	% pause()
	
	% Temporal break
	% if mode_norm(i) > 65
	% 	break
	% end
end

end
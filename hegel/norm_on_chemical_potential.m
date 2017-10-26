%% Mode with linear analog $\mu = 1$
clc; clear

% Parameters
% mu = 0.95; Omega = 1; P1 = 2;
% mu = 2.95; Omega = 8; P1 = 2;
mu = 4.95; Omega = 1; P1 = 0;

params = [mu, Omega, P1];

% Dichotomy precision
eps = 1e-9;

% Shoting from the point
xstart = -10; xend = -xstart;
xspan = [xstart, 0];

% Number of intervals
intervals = 4096;

%% Plot $u_{x, end}$ on shouting parameter $C$

cstart = 0; cstep = 0.01; cend = 1;
c = cstart:cstep:cend;

ux_end = zeros(1, length(c));
% u_end = zeros(1, length(c));

for i = 1:length(c);
	init = asympt_left(params, c(i), xstart);
	ux_end(i) = get_ux_end('f_solve', params, c(i), xspan);
	% u_end(i) = get_u_end('f_solve', params, c(i), xspan);
	
	if abs(ux_end(i)) > 100
	% if abs(u_end(i)) > 100
		break
	end
end

plot(c, ux_end); grid on
% plot(c, u_end); grid on

%% Obtaining mode manually

% mode_cspan = [0.2 0.4]; % mu = 0.95, P1 = 0, 1-st
% mode_cspan = [0.1 0.2]; % mu = 0.95, Omega = 1, P1 = 2, 1-st
% mode_cspan = [0.2 0.4]; % mu = 0.95, Omega = 8, P1 = 2, 1-st
% mode_cspan = [0.2 0.4]; % mu = 0.95, Omega = 12, P1 = 2, 1-st
% mode_cspan = [0.3 0.6]; % mu = 2.95, P1 = 0, 2-nd
% mode_cspan = [0.3 0.6]; % mu = 2.95, Omega = 8, P1 = 2, 2-nd
% mode_cspan = [0.3 0.6]; % mu = 2.95, Omega = 12, P1 = 2, 2-nd
mode_cspan = [0.3 0.6]; % mu = 4.95, Omega = 1, P1 = 0, 3-st
% mode_cspan = [0.3 0.6]; % mu = 4.95, Omega = 8, P1 = 2, 3-st
% mode_cspan = [0.3 0.6]; % mu = 4.95, Omega = 12, P1 = 2, 3-st

get_ux_end_params = @(c) get_ux_end('f_solve', params, c, xspan);
% get_u_end_params = @(c) get_u_end('f_solve', params, c, xspan);
cmode = dichotomy(get_ux_end_params, mode_cspan(1), mode_cspan(2), eps);
% cmode = dichotomy(get_u_end_params, mode_cspan(1), mode_cspan(2), eps);

[X, U] = get_symmetric_mode('f_solve', params, cmode, xspan);
% [X, U] = get_antisymmetric_mode('f_solve', params, cmode, xspan);

figure; plot_mode(params, X, U);

%% Continuation on the parameter

% mu_start = 0.95; mu_end = -4; mu_step = -0.01;
% mu_start = 2.95; mu_end = 0; mu_step = -0.01;
mu_start = 4.95; mu_end = 2.5; mu_step = -0.01;
mu = mu_start:mu_step:mu_end;

c = zeros(1, length(mu));
c(1) = cmode;

Norm = zeros(1, length(mu));
Norm(1) = get_norm(X, U);

for i = 2:length(mu)
	fprintf('%i of %i\n', i, length(mu))
	params = [mu(i) Omega P1];
	get_ux_end_params = @(c) get_ux_end('f_solve', params, c, xspan);
	% get_u_end_params = @(c) get_u_end('f_solve', params, c, xspan);
	
	c(i) = newton(get_ux_end_params, c(i - 1));
	% c(i) = newton(get_u_end_params, c(i - 1));
	[X, U] = get_symmetric_mode('f_solve', params, c(i), xspan);
	% [X, U] = get_antisymmetric_mode('f_solve', params, c(i), xspan);
	Norm(i) = get_norm(X, U);
	plot_mode(params, X, U); pause(1e-2);
end

%% 1-st mode several plots for different P1 and Omega
clc; clear

% \mu
mu_start = 0.95; mu_end = -2; mu_step = -0.01;

% Dichotomy precision
eps = 1e-9;

% Shoting from the point
xstart = -10; xend = -xstart;
xspan = [xstart, 0];

% Number of intervals
intervals = 4096;

pair_of_params = [[1 0]; [8 2]; [12 2]];
modes_cspans = [[0.2 0.4]; [0.2 0.4]; [0.2 0.4]];
mu = mu_start:mu_step:mu_end;
modes_norms_1st = zeros(3, length(mu));

for m = 1:3
	pair = pair_of_params(m, :);
	cspan = modes_cspans(m, :);
	Omega = pair(1); P1 = pair(2);
	
	params = [mu_start Omega P1];
	get_ux_end_params = @(c) get_ux_end('f_solve', params, c, xspan);
	cmode = dichotomy(get_ux_end_params, cspan(1), cspan(2), eps);
	[X, U] = get_symmetric_mode('f_solve', params, cmode, xspan);
	
	c = zeros(1, length(mu));
	c(1) = cmode;
	modes_norms_1st(m, 1) = get_norm(X, U);
	
	for i = 2:length(mu)
		fprintf('%i of %i\n', i, length(mu))
		params = [mu(i) Omega P1];
		get_ux_end_params = @(c) get_ux_end('f_solve', params, c, xspan);

		c(i) = newton(get_ux_end_params, c(i - 1));
		[X, U] = get_symmetric_mode('f_solve', params, c(i), xspan);
		modes_norms_1st(m, i) = get_norm(X, U);
	end
end

figure; hold on; grid on
plot(mu, modes_norms_1st(1, :), '-', 'LineWidth', 2, 'Color', 'black')
plot(mu, modes_norms_1st(2, :), '--', 'LineWidth', 1, 'Color', 'black')
plot(mu, modes_norms_1st(3, :), '-', 'LineWidth', 1, 'Color', 'black')

%% 2-nd mode several plots for different P1 and Omega

% \mu
mu_start = 2.95; mu_end = 0; mu_step = -0.01;

% Dichotomy precision
eps = 1e-9;

% Shoting from the point
xstart = -10; xend = -xstart;
xspan = [xstart, 0];

% Number of intervals
intervals = 4096;

pair_of_params = [[1 0]; [8 2]; [12 2]];
modes_cspans = [[0.3 0.6]; [0.3 0.6]; [0.3 0.6]];
mu = mu_start:mu_step:mu_end;
modes_norms_2nd = zeros(3, length(mu));

for m = 1:3
	pair = pair_of_params(m, :);
	cspan = modes_cspans(m, :);
	Omega = pair(1); P1 = pair(2);
	
	params = [mu_start Omega P1];
	get_u_end_params = @(c) get_u_end('f_solve', params, c, xspan);
	cmode = dichotomy(get_u_end_params, cspan(1), cspan(2), eps);
	[X, U] = get_antisymmetric_mode('f_solve', params, cmode, xspan);
	
	c = zeros(1, length(mu));
	c(1) = cmode;
	modes_norms_2nd(m, 1) = get_norm(X, U);
	
	for i = 2:length(mu)
		fprintf('%i of %i\n', i, length(mu))
		params = [mu(i) Omega P1];
		get_u_end_params = @(c) get_u_end('f_solve', params, c, xspan);

		c(i) = newton(get_u_end_params, c(i - 1));
		[X, U] = get_antisymmetric_mode('f_solve', params, c(i), xspan);
		modes_norms_2nd(m, i) = get_norm(X, U);
	end
end

plot(mu, modes_norms_2nd(1, :), '-', 'LineWidth', 2, 'Color', 'black')
plot(mu, modes_norms_2nd(2, :), '--', 'LineWidth', 1, 'Color', 'black')
plot(mu, modes_norms_2nd(3, :), '-', 'LineWidth', 1, 'Color', 'black')

%% 3-st mode several plots for different P1 and Omega

% \mu
mu_start = 4.95; mu_end = 2.5; mu_step = -0.01;

% Dichotomy precision
eps = 1e-9;

% Shoting from the point
xstart = -10; xend = -xstart;
xspan = [xstart, 0];

% Number of intervals
intervals = 4096;

pair_of_params = [[1 0]; [8 2]; [12 2]];
modes_cspans = [[0.3 0.6]; [0.3 0.6]; [0.3 0.6]];
mu = mu_start:mu_step:mu_end;
modes_norms_3rd = zeros(3, length(mu));

for m = 1:3
	pair = pair_of_params(m, :);
	cspan = modes_cspans(m, :);
	Omega = pair(1); P1 = pair(2);
	
	params = [mu_start Omega P1];
	get_ux_end_params = @(c) get_ux_end('f_solve', params, c, xspan);
	cmode = dichotomy(get_ux_end_params, cspan(1), cspan(2), eps);
	[X, U] = get_symmetric_mode('f_solve', params, cmode, xspan);
	
	c = zeros(1, length(mu));
	c(1) = cmode;
	modes_norms_3rd(m, 1) = get_norm(X, U);
	
	for i = 2:length(mu)
		fprintf('%i of %i\n', i, length(mu))
		params = [mu(i) Omega P1];
		get_ux_end_params = @(c) get_ux_end('f_solve', params, c, xspan);

		c(i) = newton(get_ux_end_params, c(i - 1));
		[X, U] = get_symmetric_mode('f_solve', params, c(i), xspan);
		modes_norms_3rd(m, i) = get_norm(X, U);
	end
end

plot(mu, modes_norms_3rd(1, :), '-', 'LineWidth', 2, 'Color', 'black')
plot(mu, modes_norms_3rd(2, :), '--', 'LineWidth', 1, 'Color', 'black')
plot(mu, modes_norms_3rd(3, :), '-', 'LineWidth', 1, 'Color', 'black')

%%
xlabel('\mu'); ylabel('N')
%% Fig. 1: Zezulins figures
% Alfimov create it for himself...?

%% Fig. 2: (N, \mu) diagrams when \Omega grows up + stability
clc; clear

% 0th mode
n_0 = 0; mu_analog_0 = 2 * n_0 + 1; mu_target_0 = -1; xstart = -8;

% 0a
Omega = 0; P1 = 0; params = [mu_target_0 Omega P1];
[mu_0a, norm_0a, stability_0a] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_0, params, xstart);

% 0b
Omega = 8; P1 = 1; params = [mu_target_0 Omega P1];
[mu_0b, norm_0b, stability_0b] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_0, params, xstart);

% 0c
Omega = 12; P1 = 1; params = [mu_target_0 Omega P1];
[mu_0c, norm_0c, stability_0c] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_0, params, xstart);

% 0d
Omega = 16; P1 = 1; params = [mu_target_0 Omega P1];
[mu_0d, norm_0d, stability_0d] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_0, params, xstart);

% 1st mode
n_1 = 1; mu_analog_1 = 2 * n_1 + 1; mu_target_1 = 0; xstart = -8;

% 1a
Omega = 0; P1 = 0; params = [mu_target_1 Omega P1];
[mu_1a, norm_1a, stability_1a] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_1, params, xstart);

% 1b
Omega = 8; P1 = 1; params = [mu_target_1 Omega P1];
[mu_1b, norm_1b, stability_1b] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_1, params, xstart);

% 1c
Omega = 12; P1 = 1; params = [mu_target_1 Omega P1];
[mu_1c, norm_1c, stability_1c] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_1, params, xstart);

% 1d
Omega = 16; P1 = 1; params = [mu_target_1 Omega P1];
[mu_1d, norm_1d, stability_1d] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_1, params, xstart);

% 2nd mode
n_2 = 2; mu_analog_2 = 2 * n_2 + 1; mu_target_2 = 3; xstart = -8;

% 2a
Omega = 0; P1 = 0; params = [mu_target_2 Omega P1];
[mu_2a, norm_2a, stability_2a] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_2, params, xstart);

% 2b
Omega = 8; P1 = 1; params = [mu_target_2 Omega P1];
[mu_2b, norm_2b, stability_2b] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_2, params, xstart);

% 2c
Omega = 12; P1 = 1; params = [mu_target_2 Omega P1];
[mu_2c, norm_2c, stability_2c] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_2, params, xstart);

% 2d
Omega = 16; P1 = 1; params = [mu_target_2 Omega P1];
[mu_2d, norm_2d, stability_2d] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_2, params, xstart);

% 3rd mode
n_3 = 3; mu_analog_3 = 2 * n_3 + 1; mu_target_3 = 5; xstart = -8;

% 3a
Omega = 0; P1 = 0; params = [mu_target_3 Omega P1];
[mu_3a, norm_3a, stability_3a] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_3, params, xstart);

% 3b
Omega = 8; P1 = 1; params = [mu_target_3 Omega P1];
[mu_3b, norm_3b, stability_3b] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_3, params, xstart);

% 3c
Omega = 12; P1 = 1; params = [mu_target_3 Omega P1];
[mu_3c, norm_3c, stability_3c] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_3, params, xstart);

% 3d
Omega = 16; P1 = 1; params = [mu_target_3 Omega P1];
[mu_3d, norm_3d, stability_3d] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_3, params, xstart);

%% Plot them all!
figure('Position', [100 100 600 270]); hold on

stability_plotter_osc(mu_0a, norm_0a, stability_0a)
stability_plotter(mu_0b, norm_0b, stability_0b)
stability_plotter(mu_0c, norm_0c, stability_0c)
stability_plotter(mu_0d, norm_0d, stability_0d)
stability_plotter_osc(mu_1a, norm_1a, stability_1a)
stability_plotter(mu_1b, norm_1b, stability_1b)
stability_plotter(mu_1c, norm_1c, stability_1c)
% stability_plotter(mu_1d, norm_1d, stability_1d)
stability_plotter_osc(mu_2a, norm_2a, stability_2a)
stability_plotter(mu_2b, norm_2b, stability_2b)
stability_plotter(mu_2c, norm_2c, stability_2c)
% stability_plotter(mu_2d, norm_2d, stability_2d)
stability_plotter_osc(mu_3a, norm_3a, stability_3a)
stability_plotter(mu_3b, norm_3b, stability_3b)
stability_plotter(mu_3c, norm_3c, stability_3c)
% stability_plotter(mu_3d, norm_3d, stability_3d)

axis([-1 7.2 0 8])

xlabel('\mu'); ylabel('N')

%% Fig. 3: (N, \mu) and many branches (antisymmetric/symmetric)
% Stability analysis for modes with linear counterpart
clc; clear

P0 = 1; P1 = 2; Omega = 8;
xspan = [-5 0];

% Now I need another approach...
% Plain of the parameters
cstart = 0.1; cstep = 0.01; cend = 50;
c_values = cstart:cstep:cend;
mu_start = 8; mu_step = -0.25; mu_end = -1;
mu_values = mu_start:mu_step:mu_end;

% Dichotomy precision
eps = 1e-9;

figure
subplot(1,2,1); hold on
subplot(1,2,2); hold on

Scan_pairs = [];

% Scaning
% -> much better (!)

figure('Position', [100 100 300 200]); hold on
plot(X, U(:, 1), 'Color', [0.6 0.6 0.6], 'LIneWidth', 2) % true
plot(x, u_approx_1, 'black') % bad approximation
plot(x, u_approx_2, 'red') % good approximation
for i = 1:length(mu_values)
	fprintf('\\mu value: %i of %i\n', i, length(mu_values));
	
	params = [mu_values(i) Omega P1];
	get_u_end_params  = @(C) get_u_end ('f_solve', params, C, xspan);
	get_ux_end_params = @(C) get_ux_end('f_solve', params, C, xspan);
	
	c_prev = c_values(1);
	
	u_end_prev  = get_u_end_params (c_prev);
	ux_end_prev = get_ux_end_params(c_prev);
	
	for j = 2:length(c_values);
		% fprintf('\tC value: %i of %i\n', j, length(c_values));
		c_next = c_values(j);
		
		u_end_next  = get_u_end_params (c_next);
		ux_end_next = get_ux_end_params(c_next);
		
		if sign(u_end_prev) ~= sign(u_end_next)
			cmode = dichotomy(get_u_end_params, c_prev, c_next, eps);
			
			if (~isempty(Scan_pairs)) && min(sqrt((Scan_pairs(:, 1) - mu_values(i)) .^ 2 + (Scan_pairs(:, 2) - cmode) .^ 2)) < 1e-2
				fprintf('\tAlready exist\n');
				c_prev = c_next;
				u_end_prev  = u_end_next;
				ux_end_prev = ux_end_next;
				continue
			end
			
			fprintf('\tNew antisymmetric mode: cmode = %g\n', cmode);
			
			[X, U] = get_antisymmetric_mode('f_solve', params, cmode, xspan);
			
			% For debug purpose
			% plot(X, U);
			% pause();
			
			mode_norm = get_norm(X, U);
			
			if mode_norm > 30
				% Too much
				fprintf('--> Norm is too big!\n');
				break
			end
			
			% Not needed?
			% fprintf('\t\tnew pair: %g, %g, %g\n', cmode, mu_values(i), modes_norm);
			
			% Not needed?
			% Scan_pairs = [Scan_pairs; [mu_values(i) modes_norm]];
			
			% Continue this solution on parameter \mu to the left and right
			Scan_pairs = [Scan_pairs; continue_mode('f_solve', mu_start, params, xspan(1), X, U, cmode, 'antisymmetric')];
			Scan_pairs = [Scan_pairs; continue_mode('f_solve', mu_end, params, xspan(1), X, U, cmode, 'antisymmetric')];
		end
		
		if sign(ux_end_prev) ~= sign(ux_end_next)
			cmode = dichotomy(get_ux_end_params, c_prev, c_next, eps);
			
			if (~isempty(Scan_pairs)) && min(sqrt((Scan_pairs(:, 1) - mu_values(i)) .^ 2 + (Scan_pairs(:, 2) - cmode) .^ 2)) < 1e-2
				fprintf('\tAlready exist\n');
				c_prev = c_next;
				u_end_prev  = u_end_next;
				ux_end_prev = ux_end_next;
				continue
			end
			
			fprintf('\tNew symmetric mode: cmode = %g\n', cmode);
			
			[X, U] = get_symmetric_mode('f_solve', params, cmode, xspan);
			
			% For debug purpose
			% plot(X, U);
			% pause();
			
			mode_norm = get_norm(X, U);
			
			if mode_norm > 30
				% Too much
				fprintf('--> Norm is too big!\n');
				break
			end
			
			% Not needed?
			% fprintf('\t\tnew pair: %g, %g, %g\n', cmode, mu_values(i), modes_norm);
			
			% Not needed?
			% Scan_pairs = [Scan_pairs; [mu_values(i) modes_norm]];
			
			% Continue this solution on parameter \mu to the left and right
			Scan_pairs = [Scan_pairs; continue_mode('f_solve', mu_start, params, xspan(1), X, U, cmode, 'symmetric')];
			Scan_pairs = [Scan_pairs; continue_mode('f_solve', mu_end, params, xspan(1), X, U, cmode, 'symmetric')];
		end
		
		c_prev = c_next;
		u_end_prev  = u_end_next;
		ux_end_prev = ux_end_next;
	end
end

% Save all values!

%%
P0 = 1; P1 = 2; Omega = 8;
xspan = [-5 0];

% Now I need another approach...
% Plain of the parameters
cstart = 0; cstep = 0.05; cend = 10;
c_values = cstart:cstep:cend;

mu = 3.4; params = [mu Omega P1];

get_u_end_params  = @(C) get_u_end ('f_solve', params, C, xspan);
get_ux_end_params = @(C) get_ux_end('f_solve', params, C, xspan);

ux_ends = zeros(1, length(c_values));

for i = 1:length(c_values)
	ux_ends(i) = get_u_end_params(c_values(i));
end

figure
plot(c_values, ux_ends)

%%

c0 = 9.6;
c1 = 10;
cmode = dichotomy(get_u_end_params, c0, c1, eps);
[X, U] = get_antisymmetric_mode('f_solve', params, cmode, xspan);
figure('Position', [100 100 300 200])
plot(X, U(:, 1), 'LineWidth', 2, 'Color', 'black')
grid on

%% Fig. 3. Stability.
clc; clear

% Commons
mu_target = -1; xstart = -6; Omega = 8; P1 = 2;

% 0th
n_0 = 0; mu_analog_0 = 2 * n_0 + 1;
params = [mu_target Omega P1];
[mu_0a, norm_0a, stability_0] = get_norm_on_chemical_potential('f_solve', mu_analog_0, params, xstart);

% 1st
n_1 = 1; mu_analog_1 = 2 * n_1 + 1;
params = [mu_target Omega P1];
[mu_1a, norm_1a, stability_1] = get_norm_on_chemical_potential('f_solve', mu_analog_1, params, xstart);

% 2nd
n_2 = 2; mu_analog_2 = 2 * n_2 + 1;
params = [mu_target Omega P1];
[mu_2, norm_2, stability_2] = get_norm_on_chemical_potential('f_solve', mu_analog_2, params, xstart);

% 3rd
n_3 = 3; mu_analog_3 = 2 * n_3 + 1;
params = [mu_target Omega P1];
[mu_3, norm_3, stability_3] = get_norm_on_chemical_potential('f_solve', mu_analog_3, params, xstart);

% 4th
n_4 = 4; mu_analog_4 = 2 * n_4 + 1;
params = [mu_target Omega P1];
[mu_4, norm_4, stability_4] = get_norm_on_chemical_potential('f_solve', mu_analog_4, params, xstart);

% 5th
n_5 = 5; mu_analog_5 = 2 * n_5 + 1;
params = [mu_target Omega P1];
[mu_5, norm_5, stability_5] = get_norm_on_chemical_potential('f_solve', mu_analog_5, params, xstart);

% 6th
n_6 = 6; mu_analog_6 = 2 * n_6 + 1;
params = [mu_target Omega P1];
[mu_6, norm_6, stability_6] = get_norm_on_chemical_potential('f_solve', mu_analog_6, params, xstart);

% 7th
n_7 = 7; mu_analog_7 = 2 * n_7 + 1;
params = [mu_target Omega P1];
[mu_7, norm_7, stability_7] = get_norm_on_chemical_potential('f_solve', mu_analog_7, params, xstart);



figure; hold on
plot(mu_0a, norm_0a)
plot(mu_1a, norm_1a)
plot(mu_2, norm_2)
plot(mu_3, norm_3)
plot(mu_4, norm_4)
plot(mu_5, norm_5)
plot(mu_6, norm_6)
plot(mu_7, norm_7)

%% Fig. 4: (N, \mu) and many branches (antisymmetric/symmetric)
% Stability analysis for modes with linear counterpart
clc; clear

P0 = 0; P1 = 1; Omega = 8;
xspan = [-8 0];

% Now I need another approach...
% Plain of the parameters
cstart = 0.1; cstep = 0.01; cend = 80;
c_values = cstart:cstep:cend;
mu_start = -1; mu_step = 0.25; mu_end = 8;
mu_values = mu_start:mu_step:mu_end;

% Dichotomy precision
eps = 1e-9;

figure
subplot(1,2,1); hold on
subplot(1,2,2); hold on

Scan_pairs = [];

% Scaning
for i = 1:length(mu_values)
	fprintf('\\mu value: %i of %i\n', i, length(mu_values));
	
	params = [mu_values(i) Omega P1];
	get_u_end_params  = @(C) get_u_end ('f_sigma_solve', params, C, xspan);
	get_ux_end_params = @(C) get_ux_end('f_sigma_solve', params, C, xspan);
	
	c_prev = c_values(1);
	
	u_end_prev  = get_u_end_params (c_prev);
	ux_end_prev = get_ux_end_params(c_prev);
	
	for j = 2:length(c_values);
		% fprintf('\tC value: %i of %i\n', j, length(c_values));
		c_next = c_values(j);
		
		u_end_next  = get_u_end_params (c_next);
		ux_end_next = get_ux_end_params(c_next);
		
		if sign(u_end_prev) ~= sign(u_end_next)
			cmode = dichotomy(get_u_end_params, c_prev, c_next, eps);
			
			if (~isempty(Scan_pairs)) && min(sqrt((Scan_pairs(:, 1) - mu_values(i)) .^ 2 + (Scan_pairs(:, 2) - cmode) .^ 2)) < 1e-2
				fprintf('\tAlready exist\n');
				c_prev = c_next;
				u_end_prev  = u_end_next;
				ux_end_prev = ux_end_next;
				continue
			end
			
			fprintf('\tNew antisymmetric mode: cmode = %g\n', cmode);
			
			[X, U] = get_antisymmetric_mode('f_sigma_solve', params, cmode, xspan);
			
			% For debug purpose
			% plot(X, U);
			% pause();
			
			mode_norm = get_norm(X, U);
			
			if mode_norm > 100
				% Too much
				fprintf('--> Norm is too big!\n');
				break
			end
			
			% Not needed?
			% fprintf('\t\tnew pair: %g, %g, %g\n', cmode, mu_values(i), modes_norm);
			
			% Not needed?
			% Scan_pairs = [Scan_pairs; [mu_values(i) modes_norm]];
			
			% Continue this solution on parameter \mu to the left and right
			Scan_pairs = [Scan_pairs; continue_mode('f_sigma_solve', mu_start, params, xspan(1), X, U, cmode, 'antisymmetric')];
			Scan_pairs = [Scan_pairs; continue_mode('f_sigma_solve', mu_end, params, xspan(1), X, U, cmode, 'antisymmetric')];
		end
		
		if sign(ux_end_prev) ~= sign(ux_end_next)
			cmode = dichotomy(get_ux_end_params, c_prev, c_next, eps);
			
			if (~isempty(Scan_pairs)) && min(sqrt((Scan_pairs(:, 1) - mu_values(i)) .^ 2 + (Scan_pairs(:, 2) - cmode) .^ 2)) < 1e-2
				fprintf('\tAlready exist\n');
				c_prev = c_next;
				u_end_prev  = u_end_next;
				ux_end_prev = ux_end_next;
				continue
			end
			
			fprintf('\tNew symmetric mode: cmode = %g\n', cmode);
			
			[X, U] = get_symmetric_mode('f_sigma_solve', params, cmode, xspan);
			
			% For debug purpose
			% plot(X, U);
			% pause();
			
			mode_norm = get_norm(X, U);
			
			if mode_norm > 100
				% Too much
				fprintf('--> Norm is too big!\n');
				break
			end
			
			% Not needed?
			% fprintf('\t\tnew pair: %g, %g, %g\n', cmode, mu_values(i), modes_norm);
			
			% Not needed?
			% Scan_pairs = [Scan_pairs; [mu_values(i) modes_norm]];
			
			% Continue this solution on parameter \mu to the left and right
			Scan_pairs = [Scan_pairs; continue_mode('f_sigma_solve', mu_start, params, xspan(1), X, U, cmode, 'symmetric')];
			Scan_pairs = [Scan_pairs; continue_mode('f_sigma_solve', mu_end, params, xspan(1), X, U, cmode, 'symmetric')];
		end
		
		c_prev = c_next;
		u_end_prev  = u_end_next;
		ux_end_prev = ux_end_next;
	end
end

% Save all values!

%%
P0 = 0; P1 = 1; Omega = 8;
xspan = [-8 0];

% Now I need another approach...
% Plain of the parameters
cstart = 0; cstep = 0.05; cend = 10;
c_values = cstart:cstep:cend;

mu = 7.5; params = [mu Omega P1];

get_u_end_params  = @(C) get_u_end ('f_sigma_solve', params, C, xspan);
get_ux_end_params = @(C) get_ux_end('f_sigma_solve', params, C, xspan);

ux_ends = zeros(1, length(c_values));

for i = 1:length(c_values)
	ux_ends(i) = get_u_end_params(c_values(i));
end

figure
plot(c_values, ux_ends)

%%
c0 = 5.4;
c1 = 5.6;
cmode = dichotomy(get_u_end_params, c0, c1, eps);
[X, U] = get_antisymmetric_mode('f_sigma_solve', params, cmode, xspan);
figure('Position', [100 100 300 200])
plot(X, U(:, 1), 'LineWidth', 2, 'Color', 'black')
grid on

%% Fig. 4. Stability.
clc; clear

% Commons
mu_target = -1; xstart = -8; Omega = 8; P1 = 1;
params = [mu_target Omega P1];

n_0 = 0; mu_analog_0 = 2 * n_0 + 1;
[mu_0, norm_0, stability_0] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_0, params, xstart);

n_1 = 1; mu_analog_1 = 2 * n_1 + 1;
[mu_1, norm_1, stability_1] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_1, params, xstart);

n_2 = 2; mu_analog_2 = 2 * n_2 + 1;
[mu_2, norm_2, stability_2] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_2, params, xstart);

n_3 = 3; mu_analog_3 = 2 * n_3 + 1;
[mu_3, norm_3, stability_3] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_3, params, xstart);

n_4 = 4; mu_analog_4 = 2 * n_4 + 1;
[mu_4, norm_4, stability_4] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_4, params, xstart);

n_5 = 5; mu_analog_5 = 2 * n_5 + 1;
[mu_5, norm_5, stability_5] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_5, params, xstart);

n_6 = 6; mu_analog_6 = 2 * n_6 + 1;
[mu_6, norm_6, stability_6] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_6, params, xstart);

n_7 = 7; mu_analog_7 = 2 * n_7 + 1;
[mu_7, norm_7, stability_7] = get_norm_on_chemical_potential('f_sigma_solve', mu_analog_7, params, xstart);

figure

stability_plotter(mu_0, norm_0, stability_0)
stability_plotter(mu_1, norm_1, stability_1)
stability_plotter(mu_2, norm_2, stability_2)
stability_plotter(mu_3, norm_3, stability_3)
stability_plotter(mu_4, norm_4, stability_4)
stability_plotter(mu_5, norm_5, stability_5)
stability_plotter(mu_6, norm_6, stability_6)
stability_plotter(mu_7, norm_7, stability_7)

%% Fig. 5. Approximation

% grey and bold: exact solution.
% black and thick thin: only one term of approximation (33).
% black and red: two terms of approximation (33).

% P0 = 0, P1 = 1.
% (1) 0th mode (n = 0), \Omega = 16, \mu = 0 (\mu_analog = 1).
% (2) 1st mdoe (n = 1), \Omega = 16, \mu = 1 (\mu_analog = 3).


clc; clear

% On a finite grid
xstart = -6; xend = 6; xstep = 0.01; xspan = [xstart xend];
x = xstart:xstep:xend;

mu = 1; Omega = 16; sigma_1 = 1; params = [mu Omega sigma_1];
n = 1; mu_analog = 2 * n + 1;

u_approx_1 = get_alfimov_approximation(params, n, x);
% u_approx_2 = get_alfimov_approximation(params, n, x);
[X_approx, U_approx] = get_another_alfimov_approximation(params, n, xspan);

% True solution
[X, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);

% -> much better (!)

figure('Position', [100 100 300 200]); hold on
plot(X, U(:, 1), 'Color', [0.6 0.6 0.6], 'LIneWidth', 3) % true
% plot(X, U(:, 1), 'Color', 'black', 'LIneWidth', 2) % true
plot(x, -u_approx_1, 'blue') % bad approximation
plot(X, U_approx(:, 1), 'red') % good approximation
% axis([-3.5 3.5 0 5])
axis([-3.5 3.5 -6 6])
xlabel('x'); ylabel('u')

%% Fig. 7. Stability
clc; clear

P0 = 0; P1 = 1; Omega = 8;
xstart = -8; xspan = [xstart 0];
mu_analog = 3;

mu_1 = 0; params_1 = [mu_1 Omega P1];
% mu_2 = 2.6; params_2 = [mu_2 Omega P1];

[X0_1, U0_1] = get_mode_with_linear_counterpart('f_solve', mu_analog, params_1, xstart);
% [X0_2, U0_2] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params_2, xstart);

[Grid_1, U_1, Norm_1] = CFDS( params_1, X0_1, U0_1(:, 1) );
% [Grid_2, U_2, Norm_2] = CFDS( params_2, X0_2, U0_2(:, 1) );

%% Spectrum
eig = get_spectrum_sigma(params_1, X0_1, U0_1, 500);
figure('Position', [100 100 400 400]); hold on
plot_spectrum(params_1, eig);

%%
figure('Position', [100 100 400 400]); hold on
plot_evolution(params_2, Grid_2, U_2)

%% Fig. 6. Eigenvalues.
clc; clear

P0 = 0; P1 = 1; Omega = 16;
xstart = -8; xspan = [xstart 0];

% 0th mode
mu_analog_0 = 1; mu_target_0 = -1.5;
params_0 = [mu_target_0 Omega P1];
[all_mu_0, all_eigs_0] = get_eigs_on_('f_sigma_solve', mu_analog_0, params_0, xstart);

% 1st mode
mu_analog_1 = 3; mu_target_1 = 0.5;
params_1 = [mu_target_1 Omega P1];
[all_mu_1, all_eigs_1] = get_eigs_on_('f_sigma_solve', mu_analog_1, params_1, xstart);

% 2nd mode
mu_analog_2 = 5; mu_target_2 = 2.5;
params_2 = [mu_target_2 Omega P1];
[all_mu_2, all_eigs_2] = get_eigs_on_('f_sigma_solve', mu_analog_2, params_2, xstart);

% 3rd mode
mu_analog_3 = 7; mu_target_3 = 4.5;
params_3 = [mu_target_3 Omega P1];
[all_mu_3, all_eigs_3] = get_eigs_on_('f_sigma_solve', mu_analog_3, params_3, xstart);

%% Plot

figure('Position', [100 100 400 300]); hold on
plot_eigs_on_chemical_potential(all_mu_0, all_eigs_0)
axis([0.5 -1.5 0 11])

figure('Position', [100 100 400 300]); hold on
plot_eigs_on_chemical_potential(all_mu_1, all_eigs_1)
axis([0.5 3 0 11])

figure('Position', [100 100 400 300]); hold on
plot_eigs_on_chemical_potential(all_mu_2, all_eigs_2)
axis([2.5 5 0 11])

figure('Position', [100 100 400 300]); hold on
plot_eigs_on_chemical_potential(all_mu_3, all_eigs_3)
axis([4.5 7 0 11])

%% (?)
P0 = 0; P1 = 1; Omega = 8;
xstart = -8; xspan = [xstart 0];
mu_analog = 7;

mu = 6.95; params = [mu Omega P1];
[X, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);

get_end = @(mex_solver_name, params, C, xspan) get_u_end('f_sigma_solve', params, C, xspan);
get_end_params = @(c) get_end('f_sigma_solve', params, c, [xstart 0]);

cstep = 0.1; cstart = 0; cend = cstart + cstep;

end_cstart = get_end_params(cstart);
end_cend = get_end_params(cend);

C = [end_cstart end_cend];

while cend < 10
	cend = cend + cstep;
	end_cend = get_end_params(cend);
	
	C = [C end_cend];
end

% nonlinear_potential = @(params, x) params(3) * cos(params(2) * x);
% spectrum = get_spectrum(params, nonlinear_potential, X, U, 300);
% figure
% plot_spectrum(params, spectrum);


%% Fig. 7. Stability - zero mean
clc; clear

P0 = 0; P1 = 1; Omega = 16;
xstart = -8; xspan = [xstart 0];
mu_analog = 7;

nonlinear_potential = @(params, x) params(3) * cos(params(2) * x);

mu = 5.5708; params = [mu Omega P1];

[X0, U0] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);
% [Grid, U, Norm] = CFDS( params, X0, U0(:, 1) );

% Spectrum
% eig = get_spectrum_sigma(params, X0, U0, 500);
eig = get_spectrum(params, nonlinear_potential, X0, U0, 300);
figure('Position', [100 100 400 400]); hold on
plot_spectrum(params, eig);

%%
figure('Position', [100 100 400 400]); hold on
plot_evolution(params, Grid, U)




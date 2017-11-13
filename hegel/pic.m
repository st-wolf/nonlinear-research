%% Fig. 1: Zezulins figures
% Alfimov create it for himself...?

%% Fig. 2: (N, \mu) diagrams when \Omega grows up + stability
% clc; clear

% 0th mode
n_0 = 0; mu_analog_0 = 2 * n_0 + 1; mu_target_0 = -2; xstart = -8;

% 0a
Omega = 0; P1 = 0; params = [mu_target_0 Omega P1];
[mu_0a, norm_0a, stability_0a] = get_norm_on_chemical_potential('f_solve', mu_analog_0, params, xstart);

% 0b
Omega = 8; P1 = 2; params = [mu_target_0 Omega P1];
[mu_0b, norm_0b, stability_0b] = get_norm_on_chemical_potential('f_solve', mu_analog_0, params, xstart);

% 0c
Omega = 12; P1 = 2; params = [mu_target_0 Omega P1];
[mu_0c, norm_0c, stability_0c] = get_norm_on_chemical_potential('f_solve', mu_analog_0, params, xstart);

% 1st mode
n_1 = 1; mu_analog_1 = 2 * n_1 + 1; mu_target_1 = -0.5; xstart = -8;

% 1a
Omega = 0; P1 = 0; params = [mu_target_1 Omega P1];
[mu_1a, norm_1a, stability_1a] = get_norm_on_chemical_potential('f_solve', mu_analog_1, params, xstart);

% 1b
Omega = 8; P1 = 2; params = [mu_target_1 Omega P1];
[mu_1b, norm_1b, stability_1b] = get_norm_on_chemical_potential('f_solve', mu_analog_1, params, xstart);

% 1c
Omega = 12; P1 = 2; params = [mu_target_1 Omega P1];
[mu_1c, norm_1c, stability_1c] = get_norm_on_chemical_potential('f_solve', mu_analog_1, params, xstart);

% 2nd mode
n_2 = 2; mu_analog_2 = 2 * n_2 + 1; mu_target_2 = 2; xstart = -8;

% 2a
Omega = 0; P1 = 0; params = [mu_target_2 Omega P1];
[mu_2a, norm_2a, stability_2a] = get_norm_on_chemical_potential('f_solve', mu_analog_2, params, xstart);

% 2b
Omega = 8; P1 = 2; params = [mu_target_2 Omega P1];
[mu_2b, norm_2b, stability_2b] = get_norm_on_chemical_potential('f_solve', mu_analog_2, params, xstart);

% 2c
Omega = 12; P1 = 2; params = [mu_target_2 Omega P1];
[mu_2c, norm_2c, stability_2c] = get_norm_on_chemical_potential('f_solve', mu_analog_2, params, xstart);

% 2d
Omega = 12; P1 = 2; params = [mu_target_2 Omega P1];
[mu_2d, norm_2d, stability_2d] = get_norm_on_chemical_potential('f_solve', mu_analog_2, params, xstart);

% Plot them all!
figure

stability_plotter(mu_0a, norm_0a, stability_0a)
stability_plotter(mu_0b, norm_0b, stability_0b)
stability_plotter(mu_0c, norm_0c, stability_0c)
stability_plotter(mu_1a, norm_1a, stability_1a)
stability_plotter(mu_1b, norm_1b, stability_1b)
stability_plotter(mu_1c, norm_1c, stability_1c)
stability_plotter(mu_2a, norm_2a, stability_2a)
stability_plotter(mu_2b, norm_2b, stability_2b)
stability_plotter(mu_2c, norm_2c, stability_2c)
stability_plotter(mu_2d, norm_2d, stability_2d)

xlabel('\mu'); ylabel('N')

%% Fig. 3: (N, \mu) and many branches (symmetric)
% Stability analysis for modes with linear counterpart
P0 = 1; P1 = 2; Omega = 8;


























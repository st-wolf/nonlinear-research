%% Preparing of different figures for the paper
% Clearing
clc; clear

% Libs
addpath ../finding_stationary_modes/

%% Potential and localized symmetric / antisymmetric modes

% U(x) = h (x^2 - x_0^2)^2
% U(0) = h x_0^4; U(\pm x_0) = 0

h = 10; x0 = 1;
U = @(x) h * (x .^ 2 - x0^2) .^ 2;

x = -1.5:0.01:1.5;

figure('Position', [100 100 300 225]); hold on;
plot(x, U(x), '-', 'LineWidth', 1, 'Color', 'k');

% Nonlinearity parameter
g = 1;

% Chemical potential
mu = 8;

% Finding the solutions using the asymptotic behaviour from the article
% Alfimov, Zezulin, Nonlinearity 20, 2007.
% There the potential U(x) a x^4 - b x^2 was considered.
a = h;
b = - 2 * (x0 ^ 2) * h;

% I use chemical potential correction to put the double-well potential U(x)
% into U(0) = 0.
mu_corrected = mu - h * (x0 ^ 2);

params = [mu_corrected, a, b, g];
xspan = [-3 0];

c_symmetric = get_symmetric_mode_parameter(params, xspan);
[X, Phi_symmetric] = get_symmetric_mode(params, c_symmetric, xspan);

zero_level = 4;
plot(X, Phi_symmetric(:, 1) + zero_level, 'LineWidth', 2, 'Color', 'k');

c_asymmetric = get_asymmetric_mode_parameter(params, xspan);
[X, Phi_asymmetric] = get_asymmetric_mode(params, c_asymmetric, xspan);

plot(X, Phi_asymmetric(:, 1) + zero_level, '--', 'LineWidth', 2, 'Color', 'k');
plot(X, zero_level * ones(1, length(X)), '-.', 'Color', 'k')

xlabel('x')

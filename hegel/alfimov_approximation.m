%% Solutions for the quantum harmonic oscillator
harmonic_oscillator_diff = @(t, y) [y(2), -(mu - (t .^ 2)) .* y(1)];

% On a finite grid
xstart = -10; xstep = 0.01; xend = 10;
x = xstart:xstep:xend;

mu = 3;

N = length(x);

algebraic_diff_operator = ...
	diag((1 / (xstep ^ 2)) * ones(1, N - 1), -1) + ...
	diag((1 / (xstep ^ 2)) * ones(1, N - 1), +1) + ...
	diag((-2 / (xstep ^ 2)) + mu - (x .^ 2));

H0_solution = (1 / sqrt(sqrt(pi))) * hermite(0, x) .* exp(-0.5 * (x .^ 2));
H1_solution = (1 / sqrt(2 * factorial(1) * sqrt(pi))) * hermite(1, x) .* exp(-0.5 * (x .^ 2));

residual = algebraic_diff_operator * (H1_solution');
plot(x, H1_solution, x, residual)

%% Compare alfimov approximation with real solution 

clc; clear

% On a finite grid
xstart = -6; xstep = 0.01; xend = 6;
x = xstart:xstep:xend;

mu = 3.1; Omega = 24; sigma_1 = 1; params = [mu Omega sigma_1];
n = 2; mu_analog = 2 * n + 1;

u_approx = get_alfimov_approximation(params, n, x);
plot(x, -u_approx)

% True solution
[X, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);

hold on
plot(X, U(:, 1))

%% Compare their spectrum

clc; clear

% On a finite grid
xstart = -6; xstep = 0.002; xend = 6;
x = xstart:xstep:xend;

mu = 0.9; Omega = 40; sigma_1 = 1; params = [mu Omega sigma_1];
n = 1; mu_analog = 2 * n + 1;

u_approx = get_alfimov_approximation(params, n, x);
[X, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);

%% Here comes da spectrum!

spec_1 = get_spectrum_sigma(params, x, u_approx, 400);
spec_2 = get_spectrum_sigma(params, X, U, 300);

%% Compare another alfimov approximation with real solution

clc; clear

% On a finite grid
xstart = -6; xend = 6; xstep = 0.01; xspan = [xstart xend];
x = xstart:xstep:xend;

mu = 0; Omega = 16; sigma_1 = 1; params = [mu Omega sigma_1];
n = 0; mu_analog = 2 * n + 1;

u_approx = get_alfimov_approximation(params, n, x);
plot(x, u_approx, 'blue');

hold on;

[X_approx, U_approx] = get_another_alfimov_approximation(params, n, xspan);
plot(X_approx, U_approx, 'black')

% True solution
[X, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);

hold on
plot(X, U(:, 1), 'red')

% -> much better (!)

%% Plot max(|u(x) - u_approx(x)|) on \sqrt{\Omega}

clc; clear

% On a finite grid
xstart = -6; xend = 6; xspan = [xstart xend];

mu = 1.1; sigma_1 = 1; n = 1; mu_analog = 2 * n + 1;
Omegas = 20:32;
Max_diffs = zeros(1, length(Omegas));

for i = 1:length(Omegas)
	fprintf('Omega: %i of %i\n', i, length(Omegas))
	
	Omega = Omegas(i);
	params = [mu Omega sigma_1];
	
	[~, U_approx] = get_another_alfimov_approximation(params, n, xspan);
	[~, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);
	
	Max_diffs(i) = max(abs(U(:, 1) - U_approx));
end

%% Plot in log-scale
figure; hold on
plot(log(sqrt(Omegas)), log(1 ./ sqrt(Omegas)))
plot(log(sqrt(Omegas)), log(Max_diffs), 'red', 'LineWidth', 2)
xlabel('log(sqrt(Omega))'); ylabel('log(max(abs(u(x) - u_{approx}(x)))')

% -> faster than 1 / \sqrt{\Omega}

%% Compare their spectrum

clc; clear

% On a finite grid
xstart = -8; xend = 8; xspan = [xstart xend];


mu = 3.5; Omega = 32; sigma_1 = 1; params = [mu Omega sigma_1];
n = 3; mu_analog = 2 * n + 1;

% [X_approx, U_approx] = get_another_alfimov_approximation(params, n, xspan);
[X, U] = get_mode_with_linear_counterpart('f_sigma_solve', mu_analog, params, xstart);

figure; hold on
plot(X, U(:, 1))
% plot(X_approx, U_approx)

%% Here comes da spectrum!

% spec_1 = get_spectrum_sigma(params, X_approx, U_approx, 500);
spec_2 = get_spectrum_sigma(params, X, U, 500);

figure;
% subplot(1, 2, 1); plot_spectrum(params, spec_1);
% subplot(1, 2, 2); plot_spectrum(params, spec_2);
plot_spectrum(params, spec_2);

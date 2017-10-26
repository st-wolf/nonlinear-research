%% Approximation of the solution for \Omega \to \infty
% integrals of 2nd and 6th power of the Hermite polynomials

% On a finite grid
xstart = -10; xstep = 0.01; xend = 10;
x = xstart:xstep:xend;

% Other parameters, equation:
% u_{xx} + (\mu - x^2) u + \sigma_1 \cos(\Omega x) u^3 = 0 
Omega = 8; sigma_1 = 1;


n = 1; mu = 1.1; mu_n = 2*n + 1; % from \mu_n to \mu. \Delta \mu = \mu_n - \mu
delta_mu = mu_n - mu;

Hn_solution = (1 / sqrt(2^n * factorial(n) * sqrt(pi))) * hermite(n, x) .* exp(-0.5 * (x .^ 2));
Hn_2nd_integral = simpson(Hn_solution .^ 2, xstep);
Hn_6th_integral = simpson(Hn_solution .^ 6, xstep);

U0 = ((2 * Hn_2nd_integral * delta_mu) / (3 * Hn_6th_integral)) ^ (1/4);

u_approx = ...
	U0 * sqrt(Omega) * Hn_solution + ...
	(U0 ^ 3) * (sigma_1 / sqrt(Omega)) * (Hn_solution .^ 3) .* cos(Omega * x);

plot(x, u_approx);

%% Compare alfimov approximation with real solution 

clc; clear

% On a finite grid
xstart = -10; xstep = 0.01; xend = 10;
x = xstart:xstep:xend;

mu = 1.1; Omega = 24; sigma_1 = 1; params = [mu Omega sigma_1];
n = 1;

u_approx = get_alfimov_approximation(params, n, x);
plot(x, -u_approx)

% True solution
[X, U] = get_mode_with_linear_counterpart(3, params, xstart);

hold on
plot(X, U(:, 1))


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

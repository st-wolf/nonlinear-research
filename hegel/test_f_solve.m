%% Fidelity test
clc; clear

omega = 4; Omega = 0.5;
params = [omega, Omega];
xspan = [0 3];
init = [0.1 0.1];

intervals = 1024;

figure; hold on; grid on

% MATLAB solver RK4
ode_params = @(t, f) ode(t, f, params);
[X, U] = RK4(ode_params, xspan, init, intervals);

plot(X, U(:, 1), 'Color', 'black', 'LineWidth', 2)

% Low level C++ solver f_solve
[X, U] = f_solve(params, xspan, init, intervals);

plot(X, U(:, 1), 'Color', 'red')

%% Speed-up test

t_start = cputime;

for i = 1:1000
	% [~, ~] = RK4(ode_params, xspan, init, intervals);
	[~, ~] = f_solve(params, xspan, init, intervals);
end

% t_elapsed_RK4 = cputime - t_start;
t_elapsed_f_solve = cputime - t_start;

% 200 times speed-up!
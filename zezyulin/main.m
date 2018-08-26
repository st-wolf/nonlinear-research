%% Playground with \alpha parameter
clc; clear

mu = -1.1; omega = sqrt(2); A = 2; params = [mu omega A];
diagram(params)
axis([-0.4 1.7 -1.4 1.4])

%% Use sewer to generate solutions plots

figure('Position', [100 100 300 200]); hold on
xlabel('x'); ylabel('u')
[X, U] = sewer(params, 1, 5); plot(X, U(:, 1), 'linewidth', 2);
axis([-4 4 0 1.6])

%% Zezuylin FIG. 1 check (e) (f) and stability
clc; clear
mu = -1.1; omega = sqrt(2); A = 2; params = [mu omega A];

figure;

[X, U] = sewer(params, 1.5, 1.5);

subplot(2, 2, 1)
hold on
plot(X, U(:, 1), 'LineWidth', 2);
plot(X, linear_potential(params, X))
plot(X, nonlinear_potential(params, X))
axis([-5 5 0 3])

eigs = get_spectrum(params, @nonlinear_potential, X, U, 64);

subplot(2, 2, 3)
plot_spectrum(params, eigs)

[X, U] = sewer(params, 0.5, 3);

subplot(2, 2, 2)
hold on
plot(X, U(:, 1), 'LineWidth', 2);
plot(X, linear_potential(params, X))
plot(X, nonlinear_potential(params, X))
axis([-5 5 0 3])

eigs = get_spectrum(params, @nonlinear_potential, X, U, 64);

subplot(2, 2, 4)
plot_spectrum(params, eigs)

%% Zezyulin FIG. 1 check (a)
clc; clear
mu_start = -4; mu_end = 1;
mu = mu_start:0.05:mu_end;

norm = zeros(3, length(mu));
mass_center = zeros(3, length(mu));
stability = zeros(3, length(mu));
c_left = zeros(3, length(mu));
c_right = zeros(3, length(mu));

c_init_left = [0.5 1.5 3];
c_init_right = [3 1.5 0.5];

omega = sqrt(2); A = 2; params = [mu(1) omega A];

for j = 1:3
    [X, U, c_left_target, c_right_target] = sewer(params,...
        c_init_left(j), c_init_right(j));
    
    c_left(j, 1) = c_left_target;
    c_right(j, 1) = c_right_target;
    
    norm(j, 1) = trapz(X, U(:, 1) .^ 2);
    mass_center(j, 1) = trapz(X, X .* (U(:, 1) .^ 2)) / norm(j, 1);

    eigs = get_spectrum(params, @nonlinear_potential, X, U, 64);
    stability(j, 1) = is_stable(eigs);
end
    
for i = 2:length(mu)
    fprintf('%i of %i', i, length(mu))
    for j = 1:3
        params = [mu(i) omega, A];
        [X, U, c_left_target, c_right_target] = sewer(params,...
            c_left(j, i - 1), c_right(j, i - 1));

        c_left(j, i) = c_left_target;
        c_right(j, i) = c_right_target;
        
        norm(j, i) = trapz(X, U(:, 1) .^ 2);
        mass_center(j, i) = trapz(X, X .* (U(:, 1) .^ 2)) / norm(j, i);

        eigs = get_spectrum(params, @nonlinear_potential, X, U, 64);
        stability(j, i) = is_stable(eigs); 
    end
end

figure('Position', [100 100 500 200]); hold on
subplot(1, 2, 1)
hold on
xlabel('\mu')
ylabel('X_c')
stability_plotter(mu, mass_center(1, :), stability(1, :))
stability_plotter(mu, mass_center(2, :), stability(2, :))
stability_plotter(mu, mass_center(3, :), stability(3, :))
axis([-4 1 -1 1])

subplot(1, 2, 2)
hold on
xlabel('\mu')
ylabel('N')
stability_plotter(mu, norm(1, :), stability(1, :))
stability_plotter(mu, norm(2, :), stability(2, :))
stability_plotter(mu, norm(3, :), stability(3, :))
axis([-4 1 0 7])

%% check if there are no modes for \omega = 0
clc; clear

mu = -1; omega = 0; A = 2; params = [mu omega A];
diagram(params)
title(sprintf('\\mu = %g, A = %g', mu, A))

% [X, U] = f_solve(params, [-10 0], asympt_left(params, 10, -10), 2 ^ 12);
% plot(X, U)

%% Using `sewer`
figure; hold on
[X, U] = sewer(params, 1.5, 1.5); plot(X, U(:, 1));
[X, U] = sewer(params, 0.5, 3); plot(X, U(:, 1));
[X, U] = sewer(params, 1, 1); plot(X, U(:, 1));

%% Spectrum
eigs = get_spectrum(params, @nonlinear_potential, X, U, 64);
plot_spectrum(params, eigs)

%% Evolution
[Grid, U, Norm] = CFDS(params, X, U(:, 1));

%% 
figure('Position', [100 100 500 400]); hold on
plot_evolution( params, Grid, U )
axis([-5 5 0 30])

%%
Up = U;
Up(abs(Up) > 1) = 1;
plot_evolution(params, Grid, Up)

%%
C = 0:0.1:30;

Ux0 = zeros(1, length(C));

for i = 1:length(C)
    init = asympt_left(params, C(i), -6);
    [X, U] = f_solve(params, [-6 0], init, 2 ^ 10);
    Ux0(i) = U(end, 2);
end

plot(C, Ux0)
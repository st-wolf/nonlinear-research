%% Playground with \alpha parameter
clc; clear

omega = 1.5; alpha = -1.34; params = [omega alpha];
diagram(params)

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
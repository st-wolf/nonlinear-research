%% Playground with \alpha parameter
clc; clear

omega = 0; alpha = 5; params = [omega alpha];
diagram(params)

%% Using `sewer`
figure; hold on
[X, U] = sewer(params, 3, 0.5); plot(X, U(:, 1));
[X, U] = sewer(params, 0.5, 3); plot(X, U(:, 1));
[X, U] = sewer(params, 1, 1);   plot(X, U(:, 1));
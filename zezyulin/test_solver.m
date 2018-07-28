%% Test MEX solver for the Kolesnikova task
clc; clear

mu = -2; omega = sqrt(2); A = 2; params = [mu omega A];
ode = @(x, u) [u(2), -(mu - 0.5 * ((omega ^ 2) * (x ^ 2))) * u(1) - (1 + A * (tanh(x) ^ 2)) * (u(1) ^ 3)];

xspan = [0 1];
u0 = [1, 1];
N = 10;

%% Run MATLAB implementation of RK4
[X, U] = RK4(ode, xspan, u0, N);

%% Run MEX implementation of RK4
[X_mex, U_mex] = f_solve(params, xspan, u0, N);

%% Comparison
max(abs(X' - X_mex))
max(abs(U - U_mex))
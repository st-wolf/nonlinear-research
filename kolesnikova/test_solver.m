%% Test MEX solver for the Kolesnikova task
clc; clear

omega = 1; alpha = 1; params = [omega alpha];
ode = @(x, u) [u(2), -(omega - (x ^ 2)) * u(1) - (1 + alpha * (x ^ 2)) * (u(1) ^ 3)];

xspan = [0 1];
u0 = [1, 1];
N = 1024;

%% Run MATLAB implementation of RK4
[X, U] = RK4(ode, xspan, u0, N);

%% Run MEX implementation of RK4
[X_mex, U_mex] = f_solve(params, xspan, u0, N);

%% Comparison
max(abs(X' - X_mex))
max(abs(U - U_mex))
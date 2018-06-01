%% Test MEX sovler for the Kucheryavy task
clc; clear

omega = 1; alpha = 0; phi = pi/2;
params = [omega, alpha, phi];
ode = @(x, u) [u(2),...
        omega * u(1) - (alpha + cos(phi) * cos(2*x) + ...
        1i * sin(phi) * sin(2 * x)) * abs(u(1))^2 * u(1)];

xspan = [0 1];
u0 = [1 + 0.5i, 1 - 2i];
N = 1024;

%% Run MATLAB implementation of RK4
[X, U] = RK4(ode, xspan, u0, N);

%% Run MEX implementation of RK4
[X_mex, U_mex] = f_solve(params, xspan, u0, N);

%% Comparison
max(abs(X' - X_mex))
max(abs(U - U_mex))
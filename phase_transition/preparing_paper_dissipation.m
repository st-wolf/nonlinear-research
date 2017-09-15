%% Libs
addpath ../finding_stationary_modes/

%% Phase in the presence of the dissipation
clc; clear

% Hamiltonian
% H = \alpha S_z^2 + \beta S_x^2 - B S_x

% System of equaiton
% \dot{S}_z = 2 \beta S_x S_y - B S_y
% \dot{S}_x = -2 \alpha S_z S_y
% \dot{S}_y = 2 (\alpha - \beta)  S_z S_x + B S_z

% Dissipation
% \alpha \to \alpha \exp(-\kappa t)
% \beta \to \beta \exp(-\kappa t)

alpha = 1; beta = 0.5; B = 0.1; R = 1;
lambda_1 = beta / alpha;
Lambda_1 = B / (alpha * R); % s = R (?)

alpha_dissip = @(t, alpha, kappa) (alpha / B) * exp(-2 * kappa * t);
beta_dissip  = @(t, beta, kappa) (beta / B) * exp(-2 * kappa * t);

% S = [S_z; S_x; S_y]
dS = @(t, S, alpha, beta, kappa) [
	2 * beta_dissip(t, beta, kappa) * S(2) * S(3) - S(3);
	-2 * alpha_dissip(t, alpha, kappa) * S(1) * S(3);
	2 * (alpha_dissip(t, alpha, kappa) - beta_dissip(t, beta, kappa)) * S(1) * S(2) + S(1)
];

% S_x^2 + S_y^2 + S_z^2 = S_0^2 = R^2
% Initial condition
Sz0 = 0; Sx0 = 0;
Sy0 = sqrt(R^2 - Sx0^2 - Sz0^2);

% Phase
% \phi = \tan^{-1} (S_y / S_x)

% (I) Phase with dissipation
% Parameter fo the dissipation
kappa = 0.05;
T_fin = 70;

% RK parameters
dS_kappa = @(t, S) dS(t, S, alpha, beta, kappa);
tspan = [0 T_fin]; S0 = [Sz0; Sx0; Sy0]; N = 2^12;
[T, S_dissip] = RK4(dS_kappa, tspan, S0, N);
phase_dissip = atan(S_dissip(:, 3) ./ S_dissip(:, 2));

figure; hold on
plot(T, phase_dissip, 'Color', 'black', 'LineWidth', 2);

% (II) Phase without dissipation
% RK parameters
dS_kappa = @(t, S) dS(t, S, alpha, beta, 0);
[T, S_barrier] = RK4(dS_kappa, tspan, S0, N);
phase = atan(S_barrier(:, 3) ./ S_barrier(:, 2));

plot(T, phase, 'Color', 'blue');

% (III) Phase in a new regime without dissipation
t_new_regime = 65; % Look at the figure
new_alpha = alpha_dissip(t_new_regime, alpha, kappa);
new_beta  = beta_dissip(t_new_regime, beta, kappa);

lambda_2 = new_beta / new_alpha;
Lambda_2 = B / (new_alpha * R);

index = find(T >= t_new_regime, 1, 'first');

% RK parameters
dS_kappa = @(t, S) dS(t, S, new_alpha, new_beta, 0);
new_S0 = S_dissip(index, :); tspan_1 = [T(index) T_fin]; tspan_2 = [T(index) 0];

[T_1, S_no_barrier_1] = RK4(dS_kappa, tspan_1, new_S0, N);
[T_2, S_no_barrier_2] = RK4(dS_kappa, tspan_2, new_S0, N);
T_2 = T_2(end:-1:1); S_no_barrier_2 = S_no_barrier_2(end:-1:1, :);

T = [T_2 T_1]; S_no_barrier = [S_no_barrier_2; S_no_barrier_1];
phase = atan(S_no_barrier(:, 3) ./ S_no_barrier(:, 2));

plot(T, phase, 'Color', 'red');

% Legend
legend( ...
	sprintf('Dissipative: \\lambda = %g, \\Lambda(0) = %g', lambda_1, Lambda_1), ...
	sprintf('Non-dissipative: \\lambda = %g, \\Lambda = %g', lambda_1, Lambda_1), ...
	sprintf('Non-dissipative: \\lambda = %g, \\Lambda = %g', lambda_2, Lambda_2))

xlabel('t'); ylabel('\phi')

%% Potential transformation
figure('Position', [100 100 325 225]); hold on

K = ellipk(sqrt(lambda_2));
x = -(1.5*K):0.01:(1.5*K);

Lambda_crit = 2 * lambda_1;

alpha_crit = B / (Lambda_crit * R);
V = @(x) alpha_crit * ( ((1/4) * (Lambda_crit ^ 2) - lambda_1 * (1 - lambda_1)) * (sn(x, lambda_1) .^ 2) ...
	- Lambda_crit * cn(x, lambda_1) ) ./ (dn(x, lambda_1) .^ 2);

plot(x, V(x), 'Color', 'black', 'LineWidth', 2);

V = @(x) alpha * ( ((1/4) * (Lambda_1 ^ 2) - lambda_1 * (1 - lambda_1)) * (sn(x, lambda_1) .^ 2) ...
	- Lambda_1 * cn(x, lambda_1) ) ./ (dn(x, lambda_1) .^ 2);

K = ellipk(sqrt(lambda_1));
x = -(1.5*K):0.01:(1.5*K);

plot(x, V(x), 'Color', 'blue')

V = @(x) new_alpha * ( ((1/4) * (Lambda_2 ^ 2) - lambda_2 * (1 - lambda_2)) * (sn(x, lambda_2) .^ 2) ...
	- Lambda_2 * cn(x, lambda_2) ) ./ (dn(x, lambda_2) .^ 2);

plot(x, V(x), 'Color', 'red')


legend(...
	sprintf('\\lambda = %g, \\Lambda = %g', lambda_1, Lambda_crit), ...
	sprintf('\\lambda = %g, \\Lambda = %g', lambda_1, Lambda_1), ...
	sprintf('\\lambda = %g, \\Lambda = %g', lambda_1, Lambda_2))

xlabel('z'); ylabel('V(z)')

%% Adiabatic approximation
omega = @(t) sqrt((2*alpha_dissip(t, alpha, kappa) - 2*beta_dissip(t, beta, kappa) + B) ...
	.* (B - 2*beta_dissip(t, beta, kappa)));

phi_0 = phase_dissip(index);
phi_prime_0 = (phase_dissip(index + 1) - phase_dissip(index-1)) / (T(index + 1) - T(index - 1));

shift = atan(-phi_prime_0 / (omega(t_new_regime) * phi_0)) - omega(t_new_regime) * t_new_regime;
amplitude = phi_0 * sqrt(omega(t_new_regime)) / cos(omega(t_new_regime) * t_new_regime + shift);

phase_adiabatic = @(t) amplitude ./ (sqrt(omega(t))) .* cos(omega(t) .* t + shift);
plot(T_1, phase_adiabatic(T_1), 'Color', 'green', 'LineWidth', 2);

%% Compare equivalent systems and their solutions
clc; clear

alpha = 1; beta = 0.5; B = 1; R = 1; kappa = 0.004; T_fin = 300; t_new_regime = 200;

alpha_dissip = @(t, alpha, kappa) alpha * exp(-kappa * t);
beta_dissip  = @(t,  beta, kappa) beta  * exp(-kappa * t);

% S = [S_z; S_x; S_y]
dS = @(t, S, alpha, beta, kappa) [
	2 * beta_dissip(t, beta, kappa) * S(2) * S(3) - B * S(3);
	-2 * alpha_dissip(t, alpha, kappa) * S(1) * S(3);
	2 * (alpha_dissip(t, alpha, kappa) - beta_dissip(t, beta, kappa)) * S(1) * S(2) + B * S(1)
];

% F0 = [S_z \phi]
dF0 = @(t, F, alpha, beta, kappa) [
	B * sqrt(1 - F(1)^2) * sin(F(2)) ...
	- 2 * beta_dissip(t, beta, kappa) * (1 - F(1)^2) * sin(F(2)) * cos(F(2));
		
	-2 * alpha_dissip(t, alpha, kappa) * F(1) ...
	+ 2 * beta_dissip(t, beta, kappa) * F(1) * cos(F(2))^2 ...
	- B * F(1) * cos(F(2)) / sqrt(1 - F(1)^2)
];

% F1 = z-linearization of F0
dF1 = @(t, F, alpha, beta, kappa) [
	B * sin(F(2)) - 2 * beta_dissip(t, beta, kappa) * sin(F(2)) * cos(F(2));
	
	-2 * alpha_dissip(t, alpha, kappa) * F(1) ...
	+ 2 * beta_dissip(t, beta, kappa) * F(1) * cos(F(2))^2 ...
	- B * F(1) * cos(F(2));
];

% F2 = \phi-linearizion of F1
dF2 = @(t, F, alpha, beta, kappa) [
	(B - 2 * beta_dissip(t, beta, kappa)) * F(2);
	(2 * beta_dissip(t, beta, kappa) - 2 * alpha_dissip(t, alpha, kappa) - B) * F(1);
];

% F3 = adiabatic approximation of F2
dF3 = @(t, F, alpha, beta, kappa) [
	F(2);
	(2 * beta_dissip(t, beta, kappa) - 2 * alpha_dissip(t, alpha, kappa) - B) ...
	* (B - 2 * beta_dissip(t, beta, kappa)) * F(1);
];

% S_x^2 + S_y^2 + S_z^2 = S_0^2 = R^2
% Initial condition
Sz0 = 0; Sx0 = 0.62;
Sy0 = sqrt(R^2 - Sx0^2 - Sz0^2);
Phi0 = atan(-Sy0 / Sx0);

% RK parameters
dS_kappa = @(t, S) dS(t, S, alpha, beta, kappa);
S_init = [Sz0; Sx0; Sy0]; tspan = [0 T_fin]; N = 2^12;
[T, S] = RK4(dS_kappa, tspan, S_init, N);
phase = atan(-S(:, 3) ./ S(:, 2));

dF0_kappa = @(t, F) dF0(t, F, alpha, beta, kappa);
F0_init = [Sz0; Phi0];
[T, F0] = RK4(dF0_kappa, tspan, F0_init, N);


index = find(T >= t_new_regime, 1, 'first');

tspan_1 = [T(index) T_fin];
% F1_init = F0(index, :);
% dF1_kappa = @(t, F) dF1(t, F, alpha, beta, kappa);
% [T1, F1] = RK4(dF1_kappa, tspan_1, F1_init, N);

% dF2_kappa = @(t, F) dF2(t, F, alpha, beta, kappa);
% [T1, F2] = RK4(dF2_kappa, tspan_1, F1_init, N);

F3_init_1 = F0(index, 2);
tmp = dF0_kappa(T(index), F0(index, :));
F3_init_2 = tmp(2);
F3_init = [F3_init_1; F3_init_2];

dF3_kappa = @(t, F) dF3(t, F, alpha, beta, kappa);
[T1, F3] = RK4(dF3_kappa, tspan_1, F3_init, N);

% Analytical solution for the adiabatic approximation
omega = @(t, alpha, beta, kappa) sqrt(...
	(B - 2 * beta_dissip(t, beta, kappa) + 2 * alpha_dissip(t, alpha, kappa)) ...
	.* (B - 2 * beta_dissip(t, beta, kappa)));

omega_kappa = @(t) omega(t, alpha, beta, kappa);
omega0 = omega_kappa(T(index));

phase_shift = atan(-F3_init(2) / (omega0 * F3_init(1)));
amplitude = F3_init(1) * sqrt(omega0) / cos(phase_shift);

solution = @(t) amplitude ./ sqrt(omega_kappa(t)) .* cos(omega_kappa(t) .* (t - T(index)) + phase_shift);

figure; hold on
% plot(T, phase, 'black', 'LineWidth', 2);
plot(T, F0(:, 2), 'red');
% plot(T1, F1(:, 2), 'green', 'LineWidth', 2);
% plot(T1, F2(:, 2), '--', 'Color', 'black');
plot(T1, F3(:, 1), '-.', 'LineWidth', 2);
plot(T1, solution(T1), 'magenta');

%% Phase portrait on the sphere
clc; clear

alpha = 1; beta = 0.5; B = 0.1;
kappa = 0.0;

df = @(t, f) [2*beta * f(2) * f(3) * exp(-2*kappa*t) - B*f(3);
	-2*alpha * f(1) * f(3) * exp(-2*kappa*t);
	2*(alpha - beta) * exp(-2*kappa*t) * f(1) * f(2) + B * f(1)];

R = 1;

figure; hold on
xlabel('S_z'); ylabel('S_x'); zlabel('S_y');

for S30 = -0.6:0.1:0.6
	S20 = sqrt(R^2 - S30^2); S10 = 0;
	f0 = [S10; S20; S30];

	[T, F] = RK4(df, [0 25], f0, 1024);
	S0 = sqrt((F(:, 1).^2) + (F(:, 2).^2) + (F(:, 3).^2));

	plot3(F(1:5:end, 1), F(1:5:end, 2), F(1:5:end, 3));
end

for S30 = -0.3:0.05:0.3
	S20 = -sqrt(R^2 - S30^2); S10 = 0;
	f0 = [S10; S20; S30];

	[T, F] = RK4(df, [0 25], f0, 1024);
	S0 = sqrt((F(:, 1).^2) + (F(:, 2).^2) + (F(:, 3).^2));

	plot3(F(1:5:end, 1), F(1:5:end, 2), F(1:5:end, 3));
end

%% Temperature (thermon period) and dissipation
clc; clear;

beta = 0.5; alpha = 1; B = 0.1; R = 1;
lambda = beta / alpha;

kappa = 0.05;

alpha_dissip = @(t) alpha * exp(-2 * kappa * t);
beta_dissip  = @(t) beta  * exp(-2 * kappa * t);
Lambda_dissip = @(t) B / (alpha_dissip(t) * R);

T = [];

t = 0; tstep = 0.01;
while true
	Lambda = Lambda_dissip(t);
	
	if Lambda > (2 * lambda)
		break
	end
	
	T = [T sqrt(lambda * (1 - lambda) + (2*lambda - 1) * Lambda / 2 - Lambda^2 / 4)];
	t = t + tstep;
end

time = 0:tstep:t;

figure('Position', [100 100 325 225]); hold on;
% plot(kappa * time, T, 'black');
plot(time, 1 ./ [T 0], 'black')
% xlabel('t'); ylabel('T / T_0')
xlabel('t'); ylabel('\tau_p / \tau_0')












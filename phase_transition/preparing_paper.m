%% Preparing of different figures for the paper
% Clearing
clc; clear

% Libs
addpath ../finding_stationary_modes/

%% Potential and localized symmetric / antisymmetric modes

% U(x) = h (x^2 - x_0^2)^2
% U(0) = h x_0^4; U(\pm x_0) = 0

h = 10; x0 = 1;
U = @(x) h * (x .^ 2 - x0^2) .^ 2;

x = -1.5:0.01:1.5;

% Nonlinearity parameter
g = 1;

% Chemical potential
% mu = 8;

xspan = [-3 0];

a = h;
b = - 2 * (x0 ^ 2) * h;

% Symmetric mode, mu = 5.7355, mu_corrected = -4.2645
% I found this chemical potential with using of my finding normalized modes
% procedure based on the asymptotic behaviour from the article Alfimov, Zezulin,
% Nonlinearity 20, 2007. There the potential U(x) a x^4 - b x^2 was considered.
mu = 5.7355;

% I use chemical potential correction to put the double-well potential U(x)
% into U(0) = 0.
mu_corrected = mu - h * (x0 ^ 2);
params = [mu_corrected, a, b, g];

c_symmetric = get_symmetric_mode_parameter(params, xspan);
[X, Phi_symmetric] = get_symmetric_mode(params, c_symmetric, xspan);

% Asymmetric mode, mu = 6.5653, mu_corrected = -3.4347
% I found this chemical potential in the same way
mu = 6.5653;

mu_corrected = mu - h * (x0 ^ 2);
params = [mu_corrected, a, b, g];

c_asymmetric = get_asymmetric_mode_parameter(params, xspan);
[X, Phi_asymmetric] = get_asymmetric_mode(params, c_asymmetric, xspan);

zero_level = 4;

% Plotting using 'plotyy'
figure('Position', [100 100 325 250]); hold on;
[h_axis, h_left, h_right] = plotyy(X, U(X), [X', X'], [Phi_symmetric(:, 1), Phi_asymmetric(:, 1)]);

set(h_left, 'Color', 'black')
set(h_right(1), 'LineWidth', 2, 'Color', 'black')
set(h_right(2), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'black')

set(h_axis(1), 'YLim', [-0.5 10.5], 'YTick', [0 5 10], 'YColor', 'black')
set(h_axis(2), 'YLim', [-2 2], 'YTick', [-1 -0.5 0 0.5 1], 'YColor', 'black')

set(h_axis(1), 'XLim', [-2 2])
set(h_axis(2), 'XLim', [-2 2])

xlabel('x')
ylabel(h_axis(1), 'U(x)')
ylabel(h_axis(2), '\Phi_{\pm}(x)')

set(h_axis(1), 'XLabel', 'x');

% Then click Edit -> Current object properties -> Undock -> Compress

%% Phase potential in elliptic coordinate with two parameters: \lambda, \Lambda
clc; clear

lambda = 0.9; Lambda = 2;
V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

K = ellipk(sqrt(lambda));
x = -(4*K):0.01:(4*K);

plot(x, V(x)); grid on

%% Rabi regime
Lambda = 100;
lambda = [0 0.25 0.5 0.9];
linewidth = [2 2 1 1];
linestyle = {'--'; '-'; '--'; '-'};

% figure('Position', [100 100 500 225]); hold on

% No legend
figure('Position', [100 100 325 225]); hold on

for i = 1:length(lambda)
	V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda(i) * (1 - lambda(i))) * (sn(x, lambda(i)) .^ 2) ...
		- Lambda * cn(x, lambda(i)) ) ./ (dn(x, lambda(i)) .^ 2);
	K = ellipk(sqrt(lambda(i)));
	x = -(2*K):0.1:(2*K);
	
	plot(x, V(x), linestyle{i}, 'Color', 'black', 'LineWidth', linewidth(i));
end

h_legend = legend(...
	  sprintf('\\lambda = %g', lambda(1)), ...
	  sprintf('\\lambda = %g', lambda(2)), ...
	  sprintf('\\lambda = %g', lambda(3)), ...
	  sprintf('\\lambda = %g', lambda(4))  ...
);

% set(h_legend, 'location', 'northeastoutside');

% Axis
plot([-5, 5], [0 0], '--', 'Color', 'k');

title(sprintf('\\Lambda = %g', Lambda));
xlabel('z');
ylabel('V_0(z)');

%% Josephson regime

% Lambda = 0.5;
lambda = [0 0.25 0.5 0.9];
linewidth = [2 2 1 1];
linestyle = {'--'; '-'; '--'; '-'};

figure('Position', [100 100 500 225]); hold on

% No legend
% figure('Position', [100 100 325 225]); hold on

for i = 1:length(lambda)
	V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda(i) * (1 - lambda(i))) * (sn(x, lambda(i)) .^ 2) ...
		- Lambda * cn(x, lambda(i)) ) ./ (dn(x, lambda(i)) .^ 2);
	K = ellipk(sqrt(lambda(i)));
	x = -(2*K):0.1:(2*K);
	
	plot(x, V(x), linestyle{i}, 'Color', 'black', 'LineWidth', linewidth(i));
end

h_legend = legend(...
	  sprintf('\\lambda = %g', lambda(1)), ...
	  sprintf('\\lambda = %g', lambda(2)), ...
	  sprintf('\\lambda = %g', lambda(3)), ...
	  sprintf('\\lambda = %g', lambda(4))  ...
);

set(h_legend, 'location', 'northeastoutside');

% Axis
plot([-2*K, 2*K], [0 0], '--', 'Color', 'k');

title(sprintf('\\Lambda = %g', Lambda));
xlabel('z');
ylabel('V_0(z)');

%% Fock regime

Lambda = 0.01;
lambda = [0 0.25 0.5 0.9];
linewidth = [2 2 1 1];
linestyle = {'--'; '-'; '--'; '-'};

figure('Position', [100 100 500 225]); hold on

% No legend
% figure('Position', [100 100 325 225]); hold on

for i = 1:length(lambda)
	V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda(i) * (1 - lambda(i))) * (sn(x, lambda(i)) .^ 2) ...
		- Lambda * cn(x, lambda(i)) ) ./ (dn(x, lambda(i)) .^ 2);
	K = ellipk(sqrt(lambda(i)));
	x = -(2*K):0.1:(2*K);
	
	plot(x, V(x), linestyle{i}, 'Color', 'black', 'LineWidth', linewidth(i));
end

h_legend = legend(...
	  sprintf('\\lambda = %g', lambda(1)), ...
	  sprintf('\\lambda = %g', lambda(2)), ...
	  sprintf('\\lambda = %g', lambda(3)), ...
	  sprintf('\\lambda = %g', lambda(4))  ...
);

set(h_legend, 'location', 'northeastoutside');

% Axis
plot([-2*K, 2*K], [0 0], '--', 'Color', 'k');

title(sprintf('\\Lambda = %g', Lambda));
xlabel('z');
ylabel('V_0(z)');

%% Phase transition diagram

figure('Position', [100 100 325 250]); hold on
axis([0 1 0 2]); xlabel('\lambda'); ylabel('\Lambda')
set(gca, 'XTick', [0 0.5 1]); set(gca, 'YTick', [0 0.5 1 2]);

Lambda = @(lambda) (1 - 16 * lambda + 16 * (lambda .^ 2) + ...
	sqrt(1 + 32 * lambda - 32 * (lambda .^ 2))) ./ (4 * (2 * lambda - 1));

% Avoid NaN at zero denominator
lambda = 0.51:0.01:1;
plot([0.5 lambda], [0 Lambda(lambda)], 'LineWidth', 2, 'Color', 'black');

lambda = 0:0.01:1;
plot(lambda, 2 * lambda, 'Color', 'black');

plot([1 1], [0 2], 'Color', 'black');
plot([0 0], [0 1], 'Color', 'black');

%% S(T) plotting
clc; clear

s = 1; % ~ Number of particles

% Potential
lambda = 0.9; Lambda = 0.1; % 1st-order transition
% lambda = 0.5; Lambda = 0.5; % 2nd-order transition

V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

% Put barrier into zero E0 = V(0)
[~, V_min] = fminsearch(V, 0); % (!)

Vb = @(x) V(x) - V_min;

x = -2:0.01:2;
plot(x, Vb(x));

V0 = Vb(0);

% Normalized
planck_const = 1;
boltzmann_const = 1;

% Energy level
E = 0.001:0.005:V0; % 1st-order transition
% E = 0.001:0.001:V0; % 2nd-order transition

T = zeros(1, length(E));
S = zeros(1, length(E)); S0 = S;
tau_p = zeros(1, length(E));

fprintf('Total number of iteration: %i\n', length(E));
for i = 1:length(E)
	fprintf('%i\n', i);
	
	[xleft, xright] = roots_near_equilibrium(@(x) -s^2 * Vb(x) + E(i));
	tau_p(i) = sqrt(2) * integral(@(x) 1 ./ sqrt(s^2 * Vb(x) - E(i)),...
		xleft, xright, 'AbsTol', 1e-6);

	T(i) = planck_const / (boltzmann_const * tau_p(i));

	S(i) = 2 * sqrt(2) * integral(@(x) sqrt(s^2 * Vb(x) - E(i)),...
		xleft, xright, 'AbsTol', 1e-6) + ...
		E(i) * tau_p(i);
	
	S0(i) = V0 * tau_p(i);
end

T_thermal = min(T):0.001:(2*max(T));
S_thermal = (planck_const * V0) ./ (boltzmann_const * T_thermal);

% Plotting
figure('Position', [100 100 650 225]);

subplot(1, 2, 1); hold on
xlabel('T'); ylabel('S_{min}(T)')

plot(T, S, 'LineWidth', 2, 'Color', 'black');
plot(T_thermal, S_thermal, '--', 'Color', 'black');

subplot(1, 2, 2); hold on
xlabel('E'); ylabel('\tau_p')

plot(E, tau_p, 'LineWidth', 2, 'Color', 'black')

%% Frequency of the small thermon oscillations
omega0 = s * sqrt(-0.5 * Lambda^2 + (2 * lambda - 1) * Lambda + 2 * lambda * (1 - lambda));



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

%% Using CFDS to calculate \gamma_{\pm \pm} numerically
clc; clear

h = 10; x0 = 1;
U = @(x) h * (x .^ 2 - x0^2) .^ 2;

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

%%

[Grid, Phi, Norm] = CFDS( params, X, Phi_symmetric(:, 1).' );









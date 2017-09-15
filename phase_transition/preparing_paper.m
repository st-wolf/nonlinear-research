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

%% gamma_{+-} with variational approximation
clc; clear

% Variational anzatz
phi_p = @(x, x0, a) (exp(-((x - x0) .^ 2) / (2 * a^2)) ...
    +exp(-((x + x0) .^ 2) / (2 * a^2)));

phi_m = @(x, x0, a) (exp(-((x - x0) .^ 2) / (2 * a^2)) ...
    -exp(-((x + x0) .^ 2) / (2 * a^2)));

get_int = @(f) integral(f, -40, 40, 'RelTol', 1e-10, 'AbsTol', 1e-6);

h = 1;
U = @(x, x0) h * ((x / x0) .^ 2 - 1) .^ 2;

a = 3;
x0 = 0:0.1:10;

gamma_pp = zeros(1, length(a));
gamma_pm = gamma_pp; gamma_mm = gamma_pp;

for i = 1:length(x0)
    phi_p_norm = sqrt( get_int(@(x) phi_p(x, x0(i), a) .^ 2) );
    phi_m_norm = sqrt( get_int(@(x) phi_m(x, x0(i), a) .^ 2) );
    
    fp = @(x) (1 / phi_p_norm) * phi_p(x, x0(i), a); % phi_p_normalized
    fm = @(x) (1 / phi_m_norm) * phi_m(x, x0(i), a); % phi_m_normalized
    
    gamma_pp(i) = get_int(@(x) (fp(x) .^ 4) );
    gamma_pm(i) = get_int(@(x) (fp(x) .^ 2) .* (fm(x) .^ 2));
    gamma_mm(i) = get_int(@(x) (fm(x) .^ 4) );
end

figure('Position', [100 100 325 225]); hold on
x = x0 / a;
plot(x, gamma_pp, x, gamma_pm, x, gamma_mm);
xlabel('x_0 / a')
ylabel('\gamma_{ij} / g')
legend('\gamma_{++}', '\gamma{+-}', '\gamma_{--}')

%% Phase potential in elliptic coordinate with two parameters: \lambda, \Lambda
clc; clear

lambda = 0.25; Lambda = 1.5;
V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

K = ellipk(sqrt(lambda));
x = -(4*K):0.01:(4*K);

figure('Position', [100 100 325 225]); hold on
plot(x, V(x), 'Color', 'black'); hold on

% x_max = inv_cn(-2 * (1 - lambda) / Lambda, lambda);
% plot([-x_max -x_max], [0 V(x_max)], 'Color', 'red');
% plot([+x_max +x_max], [0 V(x_max)], 'Color', 'red');

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

Lambda = 0.5;
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
plot([-2*K, 2*K], [0 0], '--', 'Color', 'k');

title(sprintf('\\Lambda = %g', Lambda));
xlabel('z');
ylabel('V_0(z)');

%% Phase transition diagram, combined for two barriers

figure('Position', [100 100 325 250]); hold on
axis([0 1 0 2]); xlabel('\lambda'); ylabel('\Lambda')
set(gca, 'XTick', [0 0.5 1]); set(gca, 'YTick', [0 0.5 1 2]);

Lambda = @(lambda) (1 - 16 * lambda + 16 * (lambda .^ 2) + ...
	sqrt(1 + 32 * lambda - 32 * (lambda .^ 2))) ./ (4 * (2 * lambda - 1));

% Avoid NaN at zero denominator
lambda = 0.501:0.001:1;
plot([0.5 lambda], [0 Lambda(lambda)], 'LineWidth', 2, 'Color', 'black');

lambda = 0:0.01:1;
plot(lambda, 2 * lambda, 'Color', 'black');
plot(lambda, 2 * (1 - lambda), 'Color', 'black');

plot([1 1], [0 2], 'Color', 'black');
plot([0 0], [0 1], 'Color', 'black');

%% S(T) plotting, barrier at z = 0, minima at z_{min} = \cn^{-1} (\Lambda / (2 \lambda))
clc; clear

s = 1; % ~ Number of particles
alpha = 1; % Parameter, lambda = beta / alpha
mass = 1 / (2 * alpha);

% Potential
% lambda = 0.9; Lambda = 0.1; % 1st-order transition
lambda = 0.5; Lambda = 0.5; % 2nd-order transition

% V = @(x) alpha * (s^2) * ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
%	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

 V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

% Put barrier into zero E0 = V(0)
[~, V_min] = fminsearch(V, 0); % (!)

Vb = @(x) V(x) - V_min;

% figure
% x = -100:0.01:100;
% plot(x, Vb(x));

V0 = Vb(0);

% Normalized
planck_const = 1;
boltzmann_const = 1;

% Energy level
E = linspace(0.001, V0 - 0.001, 150);

T = zeros(1, length(E));
S = zeros(1, length(E)); S0 = S;
tau_p = zeros(1, length(E));

fprintf('Total number of iteration: %i\n', length(E));
for i = 1:length(E)
	fprintf('%i\n', i);
	
	[xleft, xright] = roots_near_x0(@(x) -Vb(x) + E(i), 0);
	tau_p(i) = (1 / pi) * integral(@(x) 1 ./ sqrt(Vb(x) - E(i)),...
		xleft, xright, 'AbsTol', 1e-6);

	% T(i) = planck_const / (boltzmann_const * tau_p(i));
	T(i) = 1 / tau_p(i);

	S(i) = 2 * integral(@(x) sqrt(Vb(x) - E(i)),...
		xleft, xright, 'AbsTol', 1e-6) + ...
		pi * E(i) * tau_p(i);
	
	S0(i) = V0 * tau_p(i);
end

T_thermal = min(T):0.001:(2*max(T));
% S_thermal = (planck_const * V0) ./ (boltzmann_const * T_thermal);
S_thermal = pi * V0 ./ T_thermal;

% Plotting
figure('Position', [100 100 650 225]);

subplot(1, 2, 1); hold on
xlabel('T / T_0'); ylabel('S_{min}')

plot(T, S, 'LineWidth', 2, 'Color', 'black');
plot(T_thermal, S_thermal, '--', 'Color', 'black');

subplot(1, 2, 2); hold on
xlabel('E'); ylabel('\tau_p / \tau_0')

plot(E, tau_p, 'LineWidth', 2, 'Color', 'black')

%% Something else
P = (V0 - E) / V0;

figure('Position', [100 100 325 225]);
plot(T, sqrt(P), 'Color', 'black');
axis([0, (max(T) + max(T)*0.5), -0.1, 1.1]);

xlabel('T'); h = ylabel('P');
set(h, 'Interpreter','tex','fontsize',9) 
%% Transition temperature

% 2nd order
% Frequency of the small thermon oscillations
omega0 = 2 * alpha * s * sqrt(lambda * (1 - lambda) + 0.5 * (2 * lambda - 1) * Lambda - 0.25 * Lambda^2);
T_2nd = omega0 / (2 * pi);

subplot(1, 2, 1)
plot([T_2nd T_2nd], [min(S_thermal) max(S)]);

% 1st order
% Action on the thermon
B = 2 * s * (log((2 * sqrt(lambda) + sqrt(4 * lambda^2 - Lambda^2)) ...
	/ ((2 * sqrt(lambda) - sqrt(4 * lambda^2 - Lambda^2)))) ...
	- (Lambda / sqrt(lambda * (1 - lambda))) ...
	* atan(sqrt((1 - lambda) * (2 * lambda - Lambda) * (2 * lambda + Lambda)) / Lambda));

B_num = 2 * s * sqrt(lambda) ...
	* integral(@(x) (sqrt(1 - x .^ 2) - Delta) ./ (sqrt(1 - x .^ 2) .* (1 - lambda * (x .^ 2))), ...
	-sqrt(1 - Delta^2), +sqrt(1 - Delta^2), 'AbsTol', 1e-6);

DV = alpha * (s^2) * lambda * (1 - Lambda / (2 * lambda));
T_1st = DV / B;

% subplot(1, 2, 1)
% plot([T_1st T_1st], [-min(S_thermal) max(S)]);

%% Owerre potential
clc; clear

lambda = 0.8; Lambda = 0.2;
V_owerre = @(x) lambda * ((cn(x, lambda) - (Lambda / (2 * lambda))) .^ 2) ./ (dn(x, lambda) .^ 2);
V_me = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

[~, V_min] = fminsearch(V_me, 0); % (!)

x = -10:0.01:10;

figure; hold on;
plot(x, V_owerre(x), 'Color', 'black', 'LineWidth', 2);
plot(x, V_me(x) - V_min, 'Color', 'red');

%% Check minimum of the potential
clc; clear

lambda = 0.1; Lambda = 0.1;

V_min_exp = -(4 * lambda^2 * (1 - lambda) + Lambda^2 + Lambda^4 / (4 * lambda)) / (4 * lambda * (1 - lambda) + Lambda^2);

V = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2) - V_min_exp;

% cn_x_min = Lambda * (2 * lambda - 1) / (0.5 * Lambda^2 - 2 * lambda * (1 - lambda));
cn_x_min = Lambda / (2 * lambda);
sn_x_min = sqrt(1 - cn_x_min^2);
dn_x_min = sqrt(1 - lambda * (sn_x_min^2));

V_min = (((1/4) * Lambda^2 - lambda * (1 - lambda)) * (sn_x_min^2) - Lambda * cn_x_min) / (dn_x_min^2);

x = -3:0.01:3;
plot(x, V(x)); grid on

%% Temperature dependances, fix. lambda
% Using iterational process for the obtaining the 1-st order transition
% temparature: time-consuming computation
clc; clear
alpha = 1; s = 1;
mass = 1 / (2 * alpha);
planck_const = 1;
boltzmann_const = 1;

% 2nd order
% Frequency of the small thermon oscillations
% omega0 = @(lambda, Lambda) 2 * alpha * s * sqrt(lambda * (1 - lambda) + 0.5 * (2 * lambda - 1) * Lambda - 0.25 * Lambda^2);
% T_2nd = @(lambda, Lambda) (1 / pi) * sqrt(lambda * (1 - lambda) + 0.5 * (2 * lambda - 1) * Lambda - 0.25 * Lambda^2);
T_2nd = @(lambda, Lambda) sqrt(lambda * (1 - lambda) + 0.5 * (2 * lambda - 1) * Lambda - 0.25 * Lambda^2);


% 1st order
% Action on the thermon
B = @(lambda, Lambda) 2 * (log((2 * sqrt(lambda) + sqrt(4 * lambda^2 - Lambda^2)) ...
	/ ((2 * sqrt(lambda) - sqrt(4 * lambda^2 - Lambda^2)))) ...
	- (Lambda / sqrt(lambda * (1 - lambda))) ...
	* atan(sqrt((1 - lambda) * (2 * lambda - Lambda) * (2 * lambda + Lambda)) / Lambda));

V0 = @(lambda, Lambda) (1 - Lambda / (2 * lambda))^2;

lambda = 0.9;

% Border
Lambda_border =  (1 - 16 * lambda + 16 * (lambda .^ 2) + ...
	sqrt(1 + 32 * lambda - 32 * (lambda .^ 2))) ./ (4 * (2 * lambda - 1));

Lambda_1st = 0:0.01:Lambda_border;
t_1st = zeros(1, length(Lambda_1st));

t_1st_estimation = zeros(1, length(Lambda_1st));
for i = 1:length(Lambda_1st)
	t_1st_estimation(i) = pi * V0(lambda, Lambda_1st(i)) / B(lambda, Lambda_1st(i));
end

fprintf('1st-order transition; total number of iterations: %i\n', length(Lambda_1st));
for i = 1:length(Lambda_1st)
	fprintf('%i\n', i);
	
	V = @(x) ( ((1/4) * (Lambda_1st(i) ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda_1st(i) * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2);

	[~, V_min] = fminsearch(V, 0); % (!)

	Vb = @(x) V(x) - V_min;
	V0 = Vb(0);
	
	E = iterational_process(Vb, V0, 0.5 * V0);
	[xleft, xright] = roots_near_x0(@(x) -Vb(x) + E, 0);
	period = (1/pi) * integral(@(x) 1 ./ sqrt(Vb(x) - E),...
		xleft, xright, 'AbsTol', 1e-6);
	t_1st(i) = 1 / period;
end

Lambda_2nd = Lambda_border:0.01:(2 * lambda);
Lambda_2nd = [Lambda_2nd, 2 * lambda];
t_2nd = zeros(1, length(Lambda_2nd));

fprintf('2nd-order transition; total number of iterations: %i\n', length(Lambda_2nd));
for i = 1:length(Lambda_2nd)
	fprintf('%i\n', i);
	
	t_2nd(i) = T_2nd(lambda, Lambda_2nd(i));
end

figure('Position', [100 100 325 225]); hold on;
xlabel('\Lambda'); ylabel('T_{c}^{(1)}, T_{c}^{(2)}');
title(sprintf('\\lambda = %g', lambda))

plot(Lambda_1st, t_1st, '--', 'Color', 'k', 'LineWidth', 1);
plot(Lambda_2nd, t_2nd, 'Color', 'k', 'LineWidth', 1);

plot([Lambda_border Lambda_border], [0 t_1st(end)], 'Color', 'red');
plot(Lambda_1st, t_1st_estimation, 'Color', [0 0.5 0.2], 'LineWidth', 2)

%% Temperature dependances, fix. Lambda
clc; clear
alpha = 1; s = 1;

% 2nd order
% Frequency of the small thermon oscillations
omega0 = @(lambda, Lambda) 2 * alpha * s * sqrt(lambda * (1 - lambda) + 0.5 * (2 * lambda - 1) * Lambda - 0.25 * Lambda^2);
T_2nd = @(lambda, Lambda) omega0(lambda, Lambda) / (2 * pi);

% 1st order
% Action on the thermon
B = @(lambda, Lambda) 2 * s * (log((2 * sqrt(lambda) + sqrt(4 * lambda^2 - Lambda^2)) / ((2 * sqrt(lambda) - sqrt(4 * lambda^2 - Lambda^2)))) ...
	- (Lambda / sqrt(lambda * (1 - lambda))) * atan((sqrt(1 - lambda) * (2 * lambda - Lambda) * (2 * lambda + Lambda)) / Lambda));
V0 = @(lambda, Lambda) alpha * s^2 * (1 - Lambda / (2 * lambda))^2;
T_1st = @(lambda, Lambda) V0(lambda, Lambda) / B(lambda, Lambda);

Lambda = 0.01;
lambda = 0.5*Lambda:0.01:1;

t_1st = zeros(1, length(lambda));
t_2nd = zeros(1, length(lambda));
for i = 1:length(lambda)
	t_1st(i) = T_1st(lambda(i), Lambda);
	t_2nd(i) = T_2nd(lambda(i), Lambda);
end

figure('Position', [100 100 325 225]); hold on; grid on
xlabel('\lambda'); ylabel('T_{c}^{(1)}, T_{c}^{(2)}');
title(sprintf('\\Lambda = %g', Lambda))

plot(lambda, t_1st, '--', 'Color', 'k', 'LineWidth', 1);
plot(lambda, t_2nd, 'Color', 'k', 'LineWidth', 1);

% Border
% lambda_border = (?);
% plot([lambda_border lambda_border], [0 max(t_1st)], 'Color', 'red');

%% Testing value of action on the instanton trajectory
clc; clear

s = 1; lambda = 0.5; Lambda = 0.5;
Delta = Lambda / (2 * lambda);

B = 2 * s * (log((2 * sqrt(lambda) + sqrt(4 * lambda^2 - Lambda^2)) ...
	/ ((2 * sqrt(lambda) - sqrt(4 * lambda^2 - Lambda^2)))) ...
	- (Lambda / sqrt(lambda * (1 - lambda))) ...
	* atan(sqrt((1 - lambda) * (2 * lambda - Lambda) * (2 * lambda + Lambda)) / Lambda));

B_num = 2 * s * sqrt(lambda) ...
	* integral(@(x) (sqrt(1 - x .^ 2) - Delta) ./ (sqrt(1 - x .^ 2) .* (1 - lambda * (x .^ 2))), ...
	-sqrt(1 - Delta^2), +sqrt(1 - Delta^2), 'AbsTol', 1e-6);

%% Another phase transition: S(T), tau_p(E)
clc; clear

s = 1; alpha = 1;
mass = 1 / (2 * alpha);

planck_const = 1;
boltzmann_const = 1;

lambda = 0.9; Lambda = 0.75;
Vb = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2) - Lambda;

x_max = inv_cn(-2 * (1 - lambda) / Lambda, lambda);

V0 = Vb(x_max);

% Energy level
E = linspace(0.05, V0 - 0.001, 150);

T = zeros(1, length(E));
S = zeros(1, length(E)); S0 = S;
tau_p = zeros(1, length(E));

fprintf('Total number of iteration: %i\n', length(E));
for i = 1:length(E)
	fprintf('%i\n', i);
	
	[xleft, xright] = roots_near_x0(@(x) Vb(x) - E(i), x_max);
	tau_p(i) = sqrt(2 * mass) * integral(@(x) 1 ./ sqrt(Vb(x) - E(i)),...
		xleft, xright, 'AbsTol', 1e-6);

	T(i) = planck_const / (boltzmann_const * tau_p(i));

	S(i) = 2 * sqrt(2 * mass) * integral(@(x) sqrt(Vb(x) - E(i)),...
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

%% Check the period of the thermon

omega0 = alpha * s * sqrt((Lambda^2) / (1 - lambda) - 4 * (1 - lambda));
period = 2 * pi / omega0;

%% Temperature
clc; clear

alpha = 1; s = 1;

omega0 = @(lambda, Lambda) alpha * s * sqrt((Lambda^2) / (1 - lambda) - 4 * (1 - lambda));
T_2nd = @(lambda, Lambda) omega0(lambda, Lambda) / (2 * pi);

lambda = 0.9;
Lambda_2nd = (2 * (1 - lambda)):0.01:2;

t_2nd = zeros(1, length(Lambda_2nd));

fprintf('2nd-order transition; total number of iterations: %i\n', length(Lambda_2nd));
for i = 1:length(Lambda_2nd)
	fprintf('%i\n', i);
	
	t_2nd(i) = T_2nd(lambda, Lambda_2nd(i));
end

figure('Position', [100 100 325 225]); hold on; grid on
xlabel('\Lambda'); ylabel('T_{c}^{(2)}');
title(sprintf('\\lambda = %g', lambda))

plot(Lambda_2nd, t_2nd, 'Color', 'k', 'LineWidth', 1);

%% Chech the signs in a Taylor series
% \Lambda > 2(1 - \lambda)
clc; clear

eps = 1e-3;

lambda = 0.9;
Lambda = 1;

Vb = @(x) ( ((1/4) * (Lambda ^ 2) - lambda * (1 - lambda)) * (sn(x, lambda) .^ 2) ...
	- Lambda * cn(x, lambda) ) ./ (dn(x, lambda) .^ 2) - Lambda;
	
z_max = inv_cn(-2 * (1 - lambda) / Lambda, lambda);
		
d1f = @(x) (Vb(x + eps) - Vb(x - eps)) / eps;
d2f = @(x) (Vb(x + eps) - 2 * Vb(x) + Vb(x - eps)) / (eps^2);
d3f = @(x) (Vb(x + 2 * eps) - 2 * Vb(x + eps) + 2 * Vb(x - eps) - Vb(x - 2 * eps)) / (2 * eps^2);

%% Nice picture with potential
clc; clear

% f = @(x, y) (((x .^ 2) - 1) .^ 2) + ((y .^ 2) - 1) .^ 2;
f = @(x, y) (((x .^ 2) - 1) .^ 2) + y .^ 2 + y .^ 4;
% f = @(x, y) (((x .^ 2) - 1.2) .^ 2) - exp(0.5*y .^ 2);
% f = @(x, y) (x .^ 2 + y .^ 2 - 1) .^ 2;

x = -1.5:0.02:1.5;
y = -0.8:0.02:0.8;

[X, Y] = meshgrid(x, y);
Z = f(X, Y);
% Z(((abs(X) > 1) | (abs(Y) > 1)) & (Z > 1)) = 1;
% Z(Z > 2) = 2;
a = 1;
Z(Z > a) = a;
% Z((X + Y) < -1.3) = NaN;

surf(X, Y, Z, 'EdgeColor','none');
colormap cool
grid off













%%
clc; clear

eps = 1e-6;
omega = 4; Omega = 0.5;
params = [omega Omega];

intervals = 1024;
ode_params = @(t, f) ode(t, f, params);

psi_end = []; psi_diff_end = [];

x0 = -10;

C = 0;
cstep = 0.01;

while true
	init = [asympt(params, C, x0) asympt_diff(params, C, x0)];
	[T, Psi] = RK4(ode_params, [x0 0], init, intervals);
	
	if abs(Psi(end, 1)) > 10
		break
	else
		fprintf('psi_end = %g\n', Psi(end, 1));
		
		psi_end = [psi_end Psi(end, 1)];
		psi_diff_end = [psi_diff_end Psi(end, 2)];
		
		C = C + cstep;
	end
end

%%
clc; clear

omega = 1.3; Omega = 8;
params = [omega Omega];
xspan = [-5 0];

get_u_end_params = @(C) get_u_end(params, C, xspan);

c0 = 0.01;
cstep = 0.5;
u_end_0 = get_u_end_params(c0);

eps = 1e-6;

while true
	c1 = c0 + cstep;
	u_end_1 = get_u_end_params(c1);
		
	if sign(u_end_0) ~= sign(u_end_1)
		cmode = dichotomy(get_u_end_params, c0, c1, eps);
		break
	else
		c0 = c1;
		u_end_0 = u_end_1;
	end
end

[X, U] = get_antisymmetric_mode(params, cmode, xspan);
figure
subplot(1,2,1)
plot_mode(params, X, U)

eigs = get_spectrum(params, X, U, 128);
subplot(1,2,2);
plot_spectrum(params, eigs);

%%
[Grid, U, Norm] = CFDS(params, X, U(:, 1));

%%
N = 1000;
C = linspace(0, 6, N);

ux_end = zeros(1, N);

for i = 1:N
	ux_end(i) = get_ux_end_params(C(i));
end

plot(C, ux_end);

%%
omega = 26.17; Omega = 8;
params = [omega Omega];
xspan = [-5 0];

get_ux_end_params = @(C) get_ux_end(params, C, xspan);

c0 = 0.01;
cstep = 0.5;
ux_end_0 = get_ux_end_params(c0);

eps = 1e-6;

while true
	c1 = c0 + cstep;
	ux_end_1 = get_ux_end_params(c1);
		
	if sign(ux_end_0) ~= sign(ux_end_1)
		cmode = dichotomy(get_ux_end_params, c0, c1, eps);
		break
	else
		c0 = c1;
		ux_end_0 = ux_end_1;
	end
end

[X, U] = get_symmetric_mode(params, cmode, xspan);
figure
subplot(1,2,1)
plot_mode(params, X, U)

eigs = get_spectrum(params, X, U, 128);
subplot(1,2,2);
plot_spectrum(params, eigs);

%%

for i = -15:15
    plot([(pi + 2 * pi * i) / (2 * Omega), (pi + 2 * pi * i) / (2 * Omega)], [0, 12], 'Color', 'red', 'LineWidth', 2);
end

%%
clc; clear

omega = -26.17; Omega = 8;
params = [omega Omega];
xspan = [-1 0];

get_ux_end_params = @(C) get_ux_end(params, C, xspan);

cstep = 0.001;
C = 0:cstep:0.1;

ux_end = zeros(1, length(C));

for i = 1:length(C)
    ux_end(i) = get_ux_end_params(C(i));
end

plot(C, ux_end);
% axis([0 C(end) -1 1])

%%

eps = 1e-6;
figure
cmode = dichotomy(get_ux_end_params, 0.03, 0.08, eps);
[X, U] = get_symmetric_mode(params, cmode, xspan);
plot(X, U);

%%
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
X = S + 2*randn(size(t));

Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

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

omega = 2.9; Omega = 10;
params = [omega Omega];
xspan = [-8 0];

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

%%
eigs = get_spectrum(params, X, U, 256);
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
% Let's look on our mode (\omega = 1) while varying the parmeter \omega

omegas = omega:(-0.01):2.9;
cmodes = zeros(1, length(omegas));
cmodes(1) = cmode;

for i = 2:length(omegas);
	i
	params_shift = [omegas(i), Omega];
	get_ux_end_params_shift = @(C) get_ux_end(params_shift, C, xspan);
	cmodes(i) = newton(get_ux_end_params_shift, cmodes(i - 1));
end

plot(omegas, cmodes);

%%

for i = -15:15
    plot([pi * i / Omega, pi * i / Omega], [0, 12], 'Color', 'red', 'LineWidth', 2);
end










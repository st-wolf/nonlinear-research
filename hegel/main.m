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
% Branch with linear analog; \omega = 1.
clc; clear

omega = 0.99; Omega = 18;
params = [omega Omega];
xspan = [-10 0];

get_psix_end_params = @(C) get_psix_end(params, C, xspan);

c0 = 0.01;
cstep = 0.5;
psix_end_0 = get_psix_end_params(c0); 

eps = 1e-6;

while true
	c1 = c0 + cstep;
	psix_end_1 = get_psix_end_params(c1)
	
	if sign(psix_end_0) ~= sign(psix_end_1)
		cmode = dichotomy(get_psix_end_params, c0, c1, eps);
		break
	else
		c0 = c1;
	end
end

[X, Psi] = get_symmetric_mode(params, cmode, xspan);
plot(X, Psi)

%%
% Continue solution branch on the parameter \omega

omega_shift = omega - 0.001;
params_shift = [omega_shift, Omega];

get_psix_end_params_shift = @(C) get_psix_end(params_shift, C, xspan);
cmode_shift = newton( ...
	get_psix_end_params_shift, ...
	numerical_diff(get_psix_end_params, 1e-6), ...
	cmode);

%%
% Let's look on our mode (\omega = 1) while varying the parmeter \omega

for i = 1:100
	omega
	cmode
	
	omega_shift = omega - 0.01;
	params_shift = [omega_shift, Omega];
	
	get_psix_end_params_shift = @(C) get_psix_end(params_shift, C, xspan);
	cmode_shift = newton( ...
		get_psix_end_params_shift, ...
		numerical_diff(get_psix_end_params, 1e-6), ...
		cmode);
	
	[X, Psi] = get_symmetric_mode(params_shift, cmode_shift, xspan);
	% plot(X, Psi);
	% pause()
	
	omega = omega_shift;
	cmode = cmode_shift;
end













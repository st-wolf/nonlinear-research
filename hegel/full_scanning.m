%% Full scanning
clc; clear

omega = -5:0.1:5;
Omega = 1.5;
params = [4 Omega];
xspan = [-6 6];

cstep = 0.1;
cstart = 1.6;
cend = 2; 

% Dichotomy accuracy
eps = 1e-6;

% Target function
get_ux_end_params = @(C) get_ux_end(params, C, xspan);

N = 1000;

C = linspace(0, 5, N);
ux_end = zeros(1, N);

for i = 1:length(C)
	ux_end(i) = get_ux_end_params(C(i));
end

plot(C, ux_end);

%%
cmode = dichotomy(get_ux_end_params, 1.8, 2, eps);
[X, U] = get_mode(params, cmode, xspan);
plot(X, U);

%%
c0 = cstart;

ux0 = get_ux_end_params(c0)

for i = 1:200
	c1 = c0 + cstep;
	ux1 = get_ux_end_params(c1)
	
	if sign(ux0) ~= sign(ux1)
		cmode = dichotomy(get_ux_end_params, c0, c1, eps);
		break
	else
		c0 = c1;
		ux0 = ux1;
	end
end

[X, U] = get_mode(params, cmode, xspan);
plot(X, U);
function phi_end = get_phi_end( params, C, xspan )

xstart = xspan(1); xend = xspan(2);

% From the asymtote of the solutions
asympt_params = @(x) asympt(params, C, x);

% Derivative
eps = 1e-12;

init = [ asympt_params(xstart)
		(asympt_params(xstart + eps) - asympt_params(xstart - eps)) / eps];

ode_params = @(x, f) ode(x, f, params);
[~, Phi] = RK4(ode_params, xspan, init, 1024);

phi_end = Phi(end, 1);

end
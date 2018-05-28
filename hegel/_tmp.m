%% Fig. 3: (N, \mu) and many branches (symmetric)
% Stability analysis for modes with linear counterpart
clc; clear

P0 = 1; P1 = 2; Omega = 8;
xspan = [-8 0];

% Plain of the parameters
cstart = 0.1; cstep = 0.1; cend = 250;
c_values = cstart:cstep:cend;
mu_start = -1; mu_step = 0.1; mu_end = 8;
mu_values = mu_start:mu_step:mu_end;

Scan_pairs = [];

% Dichotomy precision
eps = 1e-9; 

% Scaning
for i = 1:length(mu_values)
	fprintf('\\mu value: %i of %i\n', i, length(mu_values));
	
	params = [mu_values(i) Omega P1];
	get_u_end_params  = @(C) get_u_end ('f_solve', params, C, xspan);
	get_ux_end_params = @(C) get_ux_end('f_solve', params, C, xspan);
	
	c_prev = c_values(1);
	
	u_end_prev  = get_u_end_params (c_prev);
	ux_end_prev = get_ux_end_params(c_prev);
	
	for j = 2:length(c_values);
		fprintf('\tC value: %i of %i\n', j, length(c_values));
		c_next = c_values(j);
		
		u_end_next  = get_u_end_params (c_next);
		ux_end_next = get_ux_end_params(c_next);
		
		if sign(u_end_prev) ~= sign(u_end_next)
			cmode = dichotomy(get_u_end_params, c_prev, c_next, eps);
			[X, U] = get_antisymmetric_mode('f_solve', params, cmode, xspan);
			
			% For debug purpose
			% plot(X, U);
			% pause();
			
			modes_norm = get_norm(X, U);
			
			fprintf('\t\tnew pair: %g, %g, %g\n', cmode, mu_values(i), modes_norm);
			Scan_pairs = [Scan_pairs; [cmode mu_values(i) modes_norm]];
			
			% if modes_norm > 30
			% 	break
			% end
		end
		
		if sign(ux_end_prev) ~= sign(ux_end_next)
			cmode = dichotomy(get_ux_end_params, c_prev, c_next, eps);
			[X, U] = get_symmetric_mode('f_solve', params, cmode, xspan);
			
			% For debug purpose
			% plot(X, U);
			% pause();
			
			modes_norm = get_norm(X, U);
			
			fprintf('\t\tnew pair: %g, %g, %g\n', cmode, mu_values(i), modes_norm);
			Scan_pairs = [Scan_pairs; [cmode mu_values(i) modes_norm]];
			
			% if modes_norm > 30
			% 	break
			% end
		end
		
		c_prev = c_next;
		u_end_prev  = u_end_next;
		ux_end_prev = ux_end_next;
	end
end

% Save all values!
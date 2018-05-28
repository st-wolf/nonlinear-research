function [ Scan_pairs ] = ...
	continue_mode( mex_solver_name, mu_target, params, xstart, X, U, cmode, symmetricity )
%
% INPUT:
%

mu_start = params(1); Omega = params(2); P1 = params(3);

if strcmp(symmetricity, 'symmetric')
	% Even solution
	get_end = @(mex_solver_name, params, C, xspan) get_ux_end(mex_solver_name, params, C, xspan);
elseif strcmp(symmetricity, 'antisymmetric')
	% Odd solution
	get_end = @(mex_solver_name, params, C, xspan) get_u_end(mex_solver_name, params, C, xspan);
end

% Continue [X, U] solution on the parameter \mu

% Step may be positive or negative depends on the mu_target and mu_start
abs_mu_step = 0.01; % why not!?

if mu_target > mu_start
	% To the right
	mu_intermediate_values = mu_start:abs_mu_step:mu_target;
else
	% To the left
	mu_intermediate_values = mu_start:(-abs_mu_step):mu_target;
end

if mu_intermediate_values(end) ~= mu_target
	mu_intermediate_values = [mu_intermediate_values mu_target];
end
	
cmode_intermediate_values = zeros(1, length(mu_intermediate_values));
cmode_intermediate_values(1) = cmode;

mode_norm = zeros(1, length(mu_intermediate_values));

% For the first mode
mode_norm(1) = get_norm(X, U);

% Plotting right here? Why not!?
subplot(1,2,1)
plot(mu_intermediate_values(1), mode_norm(1), '*', 'Color', 'red')
subplot(1,2,2)
plot(mu_intermediate_values(1), cmode_intermediate_values(1), '*', 'Color', 'red')
pause(1e-5)

Scan_pairs = [mu_intermediate_values(1) cmode_intermediate_values(1) mode_norm(1)];
for i = 2:length(mu_intermediate_values)
	% fprintf('\tContinuation: %i of %i\n', i, length(mu_intermediate_values))
	params = [mu_intermediate_values(i) Omega P1];
	
	get_end_params = @(c) get_end(mex_solver_name, params, c, [xstart 0]);
	cmode_intermediate_values(i) = newton(get_end_params, cmode_intermediate_values(i - 1));
	% cmode_intermediate_values(i) = fsolve(get_end_params, cmode_intermediate_values(i - 1));
	
	if strcmp(symmetricity, 'symmetric')
		% Even solution
		[X, U] = get_symmetric_mode(mex_solver_name, params, cmode_intermediate_values(i), [xstart 0]);
	elseif strcmp(symmetricity, 'antisymmetric')
		% Odd solution
		[X, U] = get_antisymmetric_mode(mex_solver_name, params, cmode_intermediate_values(i), [xstart 0]);
	end
	
	% Computing norm of the solution
	mode_norm(i) = get_norm(X, U);
	
	if mode_norm(i) < 1e-3
		% Zeros solution
		fprintf('Zero solution\n');
		return
	end
	
	if abs(mode_norm(i) - mode_norm(i-1)) > 1
		% Disruption
		fprintf('--> Disruption: cmode = %g, \\mu = %g\n', cmode_intermediate_values(i), mu_intermediate_values(i));
		% subplot(1,2,1)
		% plot(mu_intermediate_values(i), mode_norm(i), '*', 'Color', 'blue');
		% subplot(1,2,2)
		% plot(mu_intermediate_values(i), cmode_intermediate_values(i), '*', 'Color', 'blue');
		return
	end
	
	Scan_pairs = [Scan_pairs; [mu_intermediate_values(i) cmode_intermediate_values(i) mode_norm(i)]];
	
	% Plotting right here? Why not!?
	% plot(mu_intermediate_values(i), mode_norm(i), '.', 'Color', 'black')
	% pause(1e-5)
	
	subplot(1,2,1)
	plot(mu_intermediate_values(i), mode_norm(i), '.', 'Color', 'black')
	subplot(1,2,2)
	plot(mu_intermediate_values(i), cmode_intermediate_values(i), '.', 'Color', 'blue')
	pause(1e-5)
end

end


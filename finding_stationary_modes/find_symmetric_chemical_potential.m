function [ beta_symmetric ] = find_symmetric_chemical_potential( a, b, g )
% Find chemical potential corresponding to the nomalized symmertric
% solution of the NLS equation: \Phi_{xx} + (\beta - U(x))\Phi - g \Phi^3 = 0
% 
% INPUT:
%

% Initial estimate
beta0 = 1;

% Dependence of the symmetric mode norm on the chemical potential.
% For the sake of brevity I denote it by 'f'.
f = @(beta) compute_symmetric_mode_norm([beta a b g]);

f0 = f(beta0);

% Derivative (?)
eps = 1e-2;

% Newton method convergence
delta = 1e-7;

while true
	disp 'new iteration'
	
	df0 = (f(beta0 + eps) - f(beta0 - eps)) / eps;
	beta1 = beta0 - (f0 - 1) / df0; % (!)
	f1 = f(beta1);
	
	if abs(beta1 - beta0) < delta
		break
	else
		% New iteration
		f0 = f1;
		beta0 = beta1;
	end
	
	fprintf('beta = %g, norm = %g\n', beta0, f0);
end

beta_symmetric = beta1;

end


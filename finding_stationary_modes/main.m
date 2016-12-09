%% Finding "+" mode

% Equation
% u_{xx} + (\beta - x^2 - C \exp(-x^2))u - u^3 = 0

for beta = linspace(1, 3, 20)
	beta
	[X, Phi] = find_symmetric_mode(params, xspan);
end
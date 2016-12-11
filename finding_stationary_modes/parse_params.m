function [ beta, a, b, g ] = parse_params( params )
%
% INPUT:
%	params - 
%
% OUTPUT:
%	beta - chemical potential (energy level);
%	a, b - parameters of the double-well potential;
%	g - nonlinearity parameter.
%
% All parameters should be devided by 2 in comparison with the original GPE
% problem: i \Psi_t = -1/2 \Psi_{xx} + U(x) \Psi + g |\Psi|^2 \Psi, becouse
% of the coefficient 1/2 before the second derivative.

parsed = num2cell(params);
[beta, a, b, g] = deal(parsed{:});

end


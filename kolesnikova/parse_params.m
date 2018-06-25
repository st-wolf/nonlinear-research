function [ omega, alpha ] = parse_params( params )
%
% INPUT:
%	params - main equation parameters array (see description below)
%
% OUTPUT:
%	omega - chemical potential (energy level);
%	alpha - parameter of the nonlinear parabolic potential:
%       P(x) = 1 + \alpha x^2
%

omega = params(1); alpha = params(2);

end


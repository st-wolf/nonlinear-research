function [ mu, omega, A ] = parse_params( params )
%
% INPUT:
%	params - main equation parameters array (see description below)
%
% OUTPUT:
%	mu - chemical potential (energy level);
%   omega - parameter of the linear potential V(x) = 1/2 \omega^2  x^2
%	A - parameter of the nonlinear potential: P(x) = 1 + A \tanh^2{x}

mu = params(1); omega = params(2); A = params(3);

end


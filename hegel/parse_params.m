function [ mu, Omega ] = parse_params( params )
%
% INPUT:
%	params - 
%
% OUTPUT:
%	mu - chemical potential (energy level);
%	Omega - parameter of the nonlinear cosine potential: \cos(2 \Omega x)
%

mu = params(1); Omega = params(2);

end


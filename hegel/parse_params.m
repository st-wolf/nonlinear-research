function [ omega, Omega ] = parse_params( params )
%
% INPUT:
%	params - 
%
% OUTPUT:
%	omega - chemical potential (energy level);
%	Omega - parameter of the nonlinear cosine potential: \cos(2 \Omega x)
%

omega = params(1); Omega = params(2);

end


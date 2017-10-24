function df = ode( x, f, params )
% Stationary modes for the GPR problem
%
% INPUT:
%	

[mu, Omega] = parse_params(params);
df = [f(2), -(mu - x^2) * f(1) - (1 + 2*cos(Omega * x)) * (f(1)^3) ];

end
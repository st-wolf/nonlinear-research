function df = ode( x, f, params )
% Stationary modes for the GPR problem
%
% INPUT:
%	

[omega, Omega] = parse_params(params);
df = [f(2), -(omega - x^2) * f(1) - cos(2 * Omega * x) * (f(1)^3) ];

end
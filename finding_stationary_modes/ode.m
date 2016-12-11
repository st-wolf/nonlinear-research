function df = ode( x, f, params )
% Stationary modes for the GPE porblem
%
% INPUT:
%	

[beta, a, b, g] = parse_params(params);

df = zeros(1, 2);
df(1) = f(2);
df(2) = -(beta - potential(a, b, x)) * f(1) + g * f(1)^3;

end


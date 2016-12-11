function f = asympt( params, C, x )
% Asymtotical behaviour on the -\infty
%
% INPUT:
%	C - constant of the asymptote

[beta, a, b, ~] = parse_params(params);

f = -(C / x) * exp(...
	+(sqrt(a) / 3) * x^3 ...
	+(b / (2 * sqrt(a))) * x ...
	+(1 / (2 * sqrt(a))) * (beta^2 + (b^2) / (4*a)) / x);

end


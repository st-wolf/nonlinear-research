function f = asympt( params, C, x )
% Asymtotical behaviour on the -\infty
%
% INPUT:
%	C - constant of the asymptote

[omega, ~] = parse_params(params);

f = C * (-x)^(0.5 * (omega - 1)) * exp(-(x^2) / 2);

end
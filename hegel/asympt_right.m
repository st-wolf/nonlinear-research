function init = asympt_right( params, C, x )
% Asymtotical behaviour on the +\infty
%
% INPUT:
%	C - constant of the asymptote

[omega, ~] = parse_params(params);

func = @(x) C * (+x)^(0.5 * (omega - 1)) * exp(-(x^2) / 2);

eps = 1e-6;

f = func(x);
df = (func(x + eps) - func(x - eps)) / (2 * eps);

init = [f df];

end
function init = asympt_left( params, C, x )
% Asymtotical behaviour on the -\infty
%
% INPUT:
%	C - constant of the asymptote

[mu, omega, ~] = parse_params(params);

if omega ~= 0
    func = @(x) C * (-x)^(0.5 * (mu - 1)) * exp(-(x^2) * (omega^2) / 4);
else
    func = @(x) C * exp(sqrt(-mu) * x);
end
    
eps = 1e-12;

f = func(x);
df = (func(x + eps) - func(x - eps)) / (2 * eps);

init = [f df];

end
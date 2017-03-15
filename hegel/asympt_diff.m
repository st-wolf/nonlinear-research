function df = asympt_diff( params, C, x )

eps = 1e-6;
df = (asympt(params, C, x + eps) - asympt(params, C, x - eps)) / (2 * eps);

end


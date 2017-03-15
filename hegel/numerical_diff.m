function df = numerical_diff( f, eps )
df = @(x) (f(x + eps) - f(x - eps)) / (2 * eps);
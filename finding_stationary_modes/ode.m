function df = ode( x, f, params )

parsed = num2cell(params);
[beta, a, b] = deal(parsed{:});

df = zeros(1, 2);
df(1) = f(2);
df(2) = -(beta - a*(x .^ 4) - b*(x .^ 2)) * f(1) + f(1)^3;

end


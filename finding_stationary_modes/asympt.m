function f = asympt( params, C, x )

parsed = num2cell(params);
[beta, a, b] = deal(parsed{:});

f = -(C / x) * exp(...
	+(sqrt(a) / 3) * x^3 ...
	+(b / (2 * sqrt(a))) * x ...
	+(1 / (2 * sqrt(a))) * (beta^2 + (b^2) / (4*a)) / x);

end


function zero = dichotomy (func, a, b, eps)

fa = func(a);
while ((abs(b-a) / 2) >= eps)
	c = (a + b) / 2;
	fc = func(c);
	if sign(fc) ~= sign(fa)
		b = c;
	else
		a = c;
		fa = fc;
	end
end
zero = (a + b) / 2;

end
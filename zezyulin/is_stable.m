function stable = is_stable( eigs )

max_real = max(abs(real(eigs)));

if max_real < 1e-2
	stable = true;
else
	stable = false;
end

end


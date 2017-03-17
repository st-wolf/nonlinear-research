function q = simpson( f, dx )
q = (dx / 3) * sum(f(1:2:end-2) + 4*f(2:2:end-1) + f(3:2:end));


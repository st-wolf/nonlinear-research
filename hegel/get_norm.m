function Norm = get_norm(X, U)
dx = X(2) - X(1);
Norm = simpson(U(:, 1) .^ 2, dx);
end


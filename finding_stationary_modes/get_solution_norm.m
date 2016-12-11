function [ phi_norm ] = get_solution_norm( X, Phi )
% Compute norm of the solutions using the Simpson's rule
%

dx = X(2) - X(1);
phi_norm = sqrt(simpson(Phi(:, 1) .^ 2, dx));

end


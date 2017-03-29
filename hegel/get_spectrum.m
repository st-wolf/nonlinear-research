function eigenvalues = get_spectrum( params, X, U, N )

% TODO Symmetrize with respect to the horizontal axis

% Stability of localized nonlinear modes
% Here we use Jianke Yang method

% INPUT:
% :params: parameter of NLS with nonlinear potential (see info.txt)
% :X: grid
% :U: values on grid
% :N: number of harmonics in fourie series

% OUTPUT:
% :eigenvalues: eigen values of differential operator, that correspon to
%	the evolution of perturbations

period = abs(X(end) - X(1));
xstep = X(2) - X(1);

% Wavenumber
k = 2*pi / period;

G{1} = params(1) - (X .^ 2) + f_potential(params, X) .* (U(:, 1).^2);
G{2} = params(1) - (X .^ 2) + 3*f_potential(params, X) .* (U(:, 1).^2);

for n = -N:N
	Gf{1}(n+N+1) = (1 / period) * simpson(G{1} .* exp(-1i * k*n*X), xstep);
	Gf{2}(n+N+1) = (1 / period) * simpson(G{2} .* exp(-1i * k*n*X), xstep);
end

% Big Gf{:}(end) values indicate bad convergence of the Fourie series
fprintf('Spectrum localization: Gf(end) = %g, %g\n', Gf{1}(end), Gf{2}(end));

D = -k^2 * (diag(-N:N) .^ 2);

Gt{1} = toeplitz([Gf{1}(N+1:2*N+1) zeros(1, N)],...
					  [Gf{1}(N+1:-1:1)  zeros(1, N)]);
Gt{2} = toeplitz([Gf{2}(N+1:2*N+1) zeros(1 ,N)],...
					  [Gf{2}(N+1:-1:1)  zeros(1 ,N)]);

M = (D + Gt{2}) * (D + Gt{1});

eigenvalues = sqrt(-eig(M));
eigenvalues = [eigenvalues; -eigenvalues];

end

% ------------------------------------------------------------------------
function potential = f_potential(params, x)
potential = cos(2 * params(2) * x);
end

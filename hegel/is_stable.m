function [ stable ] = is_stable( params, nonlinear_potential, X, U, N )
% Use @get_spectrum to determine the stability of the mode
%
% INPUT:
% :params: parameter of NLS with nonlinear potential (see info.txt)
% :nonlinear_potential: function of the form @nonlinear_potential(params, x)
% :X: grid
% :U: values on grid
% :N: number of harmonics in fourie series

% OUTPUT:
% :stable: true of false
%

eigenvalues = get_spectrum(params, nonlinear_potential, X, U, N);
max_real = max(abs(real(eigenvalues)));

% For debug purpoces
% plot_spectrum(params, eigenvalues);

% Compareson with zero
if max_real < 1e-3
	stable = true;
else
	stable = false;

end


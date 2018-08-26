function plot_spectrum( params, eigenvalues )

% INPUT:
% :params: parameter of NLS with nonlinear potential (see info.txt)
% :eigenvalues: eigen values of differential operator, that correspon to
%	the evolution of perturbations

hold on

% Continuous spectrum
% line([0 0], [ params(1)  100], 'Color', 'r', 'LineWidth', 2)
% line([0 0], [-params(1) -100], 'Color', 'r', 'LineWidth', 2)

% Discrete eigenvalues
plot(eigenvalues, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 5)

% Axis boundaries
xmax = max(1, 1.2 * max(real(eigenvalues)));
ymax = 1.2 * max(imag(eigenvalues(real(eigenvalues) > 1e-5)));

if ~isempty(ymax)
	ymax = max(15, ymax);
else
	ymax = 15;
end

axis([-xmax xmax -ymax ymax])
% axis([-xmax xmax -3 3])

% title(sprintf('Parameters: \\mu = %g, \\omega = %g, A = %g', params(1), params(2), params(3)))
grid on

end


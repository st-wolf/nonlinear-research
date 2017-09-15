function plot_mode( params, X, U )

% INPUT:
% :params: parameter of NLS with nonlinear potential (see info.txt)
% :X: grid
% :U: values on grid

% figure
hold on

umax = max(U(:, 1));
umin = min(U(:, 1));

%line([-2*pi:pi:2*pi; -2*pi:pi:2*pi],...
%	  [repmat(-2*abs(umin), 1, 5); repmat(2*abs(umax), 1, 5)],...
%	  'Color', 'r', 'LineWidth', 1)

plot(X, U(:, 1), '-k', 'LineWidth', 2)

xlim([-2*pi, 2*pi])
ylim([1.2 * umin, 1.2 * umax])

title(sprintf('Parameters: omega = %g, Omega = %g', params(1), params(2)))
grid on
hold off
end


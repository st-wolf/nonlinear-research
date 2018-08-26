function plot_evolution( params, Grid, U )

xrange = [-5 5];
X = Grid{1}(1, :);
T = Grid{2}(:, 1);

sel = ((X > xrange(1)) & (X < xrange(2)));

image([xrange(1) xrange(2)], [T(1) T(end)], abs(U(:, sel)), 'CDataMapping', 'scaled')
set(gca,'YDir','normal');
xlabel('X')
ylabel('T')
% title_str = sprintf('Parameters: omega = %g, alpha = %g', params(1), params(2));
% title(title_str)

end


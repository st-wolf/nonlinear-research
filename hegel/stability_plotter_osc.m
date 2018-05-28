function stability_plotter_osc( X, Y, stability )
%
% INPUT:
%	:X:

hold on;

for i = 2:length(X)
	if stability(i - 1)
		plot([X(i - 1) X(i)], [Y(i - 1) Y(i)], 'LineWidth', 4, 'Color', [0.6 0.6 0.6])
	else
		plot([X(i - 1) X(i)], [Y(i - 1) Y(i)], 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
	end
end

end



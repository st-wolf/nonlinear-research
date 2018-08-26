function stability_plotter( X, Y, stability )
%
% INPUT:
%	:X:

hold on;

for i = 2:length(X)
	if stability(i - 1)
		plot([X(i - 1) X(i)], [Y(i - 1) Y(i)], 'LineWidth', 2, 'Color', 'black')
	else
		plot([X(i - 1) X(i)], [Y(i - 1) Y(i)], 'LineWidth', 1, 'Color', 'red')
	end
end

end


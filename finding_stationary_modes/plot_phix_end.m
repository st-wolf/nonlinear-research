function plot_phix_end( params, xspan )

C = linspace(4, 6, 100);
Phix_end = zeros(1, length(C));

for i = 1:length(C)
	Phix_end(i) = get_phix_end(params, C(i), xspan);
end

plot(C, Phix_end)

end
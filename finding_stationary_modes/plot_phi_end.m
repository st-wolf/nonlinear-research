function plot_phi_end( params, xspan )

C = linspace(0, 2, 100);
Phi_end = zeros(1, length(C));

for i = 1:length(C)
	Phi_end(i) = get_phi_end(params, C(i), xspan);
end

plot(C, Phi_end)

end

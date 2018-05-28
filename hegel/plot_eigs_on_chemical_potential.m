function plot_eigs_on_chemical_potential( mu, eigs )

hold on; set(gca, 'XDir','reverse')

for i = 1:length(mu)
	if abs(real(eigs(i))) < 0.01
		plot(mu(i), abs(eigs(i)), '.', 'Color', 'black')
	else
		plot(mu(i), abs(eigs(i)), '.', 'Color', 'red')
	end
end

end


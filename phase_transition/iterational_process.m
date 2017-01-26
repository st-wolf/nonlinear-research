function [ E ] = iterational_process( Vb, V0, mass, E0 )

eps = 1e-6;

while true
	[xleft, xright] = roots_near_equilibrium(@(x) -Vb(x) + E0);
	period = sqrt(2 * mass) * integral(@(x) 1 ./ sqrt(Vb(x) - E0),...
		xleft, xright, 'AbsTol', 1e-6);

	E = V0 - (2 * sqrt(2 * mass) / period) * integral(@(x) sqrt(Vb(x) - E0),...
		xleft, xright, 'AbsTol', 1e-6);
	
	if abs(E - E0) < eps
		return
	else
		E0 = E;
	end
end

end


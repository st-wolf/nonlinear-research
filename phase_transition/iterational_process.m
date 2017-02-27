function [ E ] = iterational_process( Vb, V0, E0 )

eps = 1e-6;

while true
	[xleft, xright] = roots_near_x0(@(x) -Vb(x) + E0, 0);
	% [xleft, xright] = roots_near_equilibrium(@(x) -Vb(x) + E0);
	period = integral(@(x) 1 ./ sqrt(Vb(x) - E0),...
		xleft, xright, 'AbsTol', 1e-6);

	E = V0 - (2 / period) * integral(@(x) sqrt(Vb(x) - E0),...
		xleft, xright, 'AbsTol', 1e-6);
	
	if abs(E - E0) < eps
		return
	else
		E0 = E;
	end
end

end


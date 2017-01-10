function [Grid, U, Norm] = CFDS( params, X0, U0 )
% Solving the equation: i u_t = -u_{xx} - (\alpha + \cos 2x) |u|^2 u

% Accuracy of interation process
eps = 1e-12;

[~, a, b, gN] = parse_params(params);

% Zero linear potential
% V = @(x) 0;

% Double-well linear potential
V = @(x) potential(a, b, x);

% Nonlinear cosine potential
% g = @(x) alpha + cos(2 * x);

% Repulsive nonlinearity
g = @(x) -gN;

xstep = 0.03;
xspan = [X0(1), X0(end)];
xgrid = xspan(1):xstep:xspan(2);
Nx = length(xgrid);

tstep = 0.01;
tspan = [0 100];
tgrid = tspan(1):tstep:tspan(2);
Nt = length(tgrid);

[Grid{1}, Grid{2}] = meshgrid(xgrid, tgrid);
U = zeros(Nt, Nx);

% Initial condition
U(1, :) = resample(X0, U0, xgrid);

% Perturbation (?)

U(1, 1) = 0;
U(1, end) = 0;

figure
plot(U(1, :))

% Norm of the solution
Norm = zeros(1, Nt);

Norm(1) = simpson(abs(U(1, :)) .^ 2, xstep);

% Prepare tridiagonal sparse linear operator: Lv = F(v, v')
	
% -1 diagonal
Ljm1 = 0.5 / (xstep ^ 2) * ones(1, Nx - 2);
Ljm1 = [Ljm1 0];
	
% +1 diagonal
Ljp1 = 0.5 / (xstep ^ 2) * ones(1, Nx - 2);
Ljp1 = [0 Ljp1];

% Dissipation (?)
% Dissip = 1i * dissipation(xgrid);

for i = 2:Nt
	% 0 diagonal
	Lj = 1i / tstep - 1 / (xstep ^ 2) - 0.5 * (V(xgrid) - 0.5 * g(xgrid) .* (abs(U(i - 1, :)) .^ 2));
	
	% Without dissipation
	L = diag(Lj) + diag(Ljm1, -1) + diag(Ljp1, +1);
	
	% With dissipation
	% L = diag(Lj) + diag(Dissip) + diag(Ljm1, -1) + diag(Ljp1, +1);
	
	% Prepare nonlinear right part: Lv = F(v, v')
	% Independent from iterations
	
	F0 = -0.5 / (xstep ^ 2) * (U(i - 1, 1:end - 2) - 2 * U(i - 1, 2:end - 1) + U(i - 1, 3:end)) + 1i / tstep * U(i - 1, 2:end - 1) + ...
		+ 0.5 * U(i - 1, 2:end - 1) .* (V(xgrid(2:end - 1)) - 0.5 * g(xgrid(2:end - 1)) .* abs(U(i - 1, 2:end - 1)) .^ 2);
	
	% Iterations
	U(i, :) = U(i - 1, :);
	
	while true
		F = F0 - 0.25 * g(xgrid(2:end - 1)) .* (abs(U(i, 2:end - 1)) .^ 2) .* (U(i, 2:end - 1) + U(i - 1, 2:end - 1));
		
		% Zero boundary condition
		F = [0; F.'; 0];
		
		L = sparse(L);
		Unext = L \ F;
		Unext = Unext.';
		
		if max(abs(Unext)) > 1e10
			fprintf('Collapsed!!!\n')
			return
		end
		
		if ( max(abs(Unext - U(i, :))) - (1e-9) * max(abs(U(i, :))) ) < eps
			U(i, :) = Unext;
			break
		else
			U(i, :) = Unext;
		end
		
	end
	
	plot(abs(U(i, :))); drawnow
	
	pause(1e-2);
	
	Norm(i) = simpson(abs(U(i, :)) .^ 2, xstep);
	
	% Logging
	fprintf('iter = %i, max = %g, norm = %g\n', i, max(abs(U(i, :))), Norm(i));
	
end

end
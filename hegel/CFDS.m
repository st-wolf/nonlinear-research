function [Grid, U, Norm] = CFDS( params, X0, U0 )
% Solving the equation:
% i \Psi_t = -\Psi_{xx} + x^2 \Psi + \cos(2 \Omega x) |\Psi|^2 \Psi

% Accuracy of interation process
eps = 1e-12;

% Linear potential
V = @(x) x .^ 2;

% Nonlinear cosine potential parameter
Omega = params(2);
P1 = params(3);

% Nonlinear cosine potential
% g = @(x) P1 * cos(Omega * x);
g = @(x) 1 + P1 * cos(Omega * x);

xstep = 0.02;
xspan = [-6 6];
xgrid = xspan(1):xstep:xspan(2);
Nx = length(xgrid);

tstep = 0.005;
tspan = [0 12];
tgrid = tspan(1):tstep:tspan(2);
Nt = length(tgrid);

[Grid{1}, Grid{2}] = meshgrid(xgrid, tgrid);
U = zeros(Nt, Nx);

% Initial condition
U(1, :) = resample(X0, U0, xgrid);

% Perturbation
a = 0.01;
f = 10;
U(1, :) = U(1, :) .* (1 + a * (1 + sin(f * xgrid)));

U(1, 1) = 0;
U(1, end) = 0;

% figure
% plot(U(1, :))

% Norm of the solution
Norm = zeros(1, Nt);
Norm(1) = simpson(abs(U(1, :)) .^ 2, xstep);

% TODO: add energy calculation

% Prepare tridiagonal sparse linear operator: Lv = F(v, v')
	
% -1 diagonal
Ljm1 = 0.5 / (xstep ^ 2) * ones(1, Nx - 2);
Ljm1 = [Ljm1 0];
	
% +1 diagonal
Ljp1 = 0.5 / (xstep ^ 2) * ones(1, Nx - 2);
Ljp1 = [0 Ljp1];


Dissip = 1i * dissipation(xgrid);

for i = 2:Nt
	% 0 diagonal
	Lj = 1i / tstep - 1 / (xstep ^ 2) - 0.5 * (V(xgrid) - 0.5 * g(xgrid) .* (abs(U(i - 1, :)) .^ 2));
	L = diag(Lj) + diag(Dissip) + diag(Ljm1, -1) + diag(Ljp1, +1);
	% L = diag(Lj) + diag(Ljm1, -1) + diag(Ljp1, +1);
	
	
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
    pause(1e-1)
	
    Norm(i) = simpson(abs(U(i, :)) .^ 2, xstep);
	% Logging
	fprintf('iter = %i of %i, max = %g, norm = %g\n', i, Nt, max(abs(U(i, :))), Norm(i));
	
end

end


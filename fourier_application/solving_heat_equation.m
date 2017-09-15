%%
clc
clear

% Parameter of the heat transition
c = 1;

xstart = -100; xend = 100; xstep = 0.01;
xgrid = xstart:xstep:xend;

tstep = 0.005;
tgrid = 0:0.01:10;

% Initial condition
phi = @(x) exp(-x .^ 2);

U = zeros(length(tgrid), length(xgrid));
U(1, :) = phi(xgrid);

% Harmonics multiplicator
N = length(xgrid);
k = (2 * pi / (xend - xstart)) * [0:(N-1)/2, -(N-1)/2:(-1)];

% For every temporal layer
for i = 2:length(tgrid)
	initials = U(i-1, :);
	initials_hat = fft(initials);
	solution_hat = initials_hat .* exp(-c * (k .^ 2) .* tstep);
	
	U(i, :) = ifft(solution_hat);
	
	% Animation
	plot(xgrid, U(i, :), 'LineWidth', 2, 'Color', 'black');
	axis([-10 10 0 1]);
	grid on;
	pause(1e-1);
end
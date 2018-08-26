function diagram( params )
% TODO: description
%
% INPUT:
%	params - main equation parameters array (see `parse_params`)

% INTERNAL PARAMETERS:
xstart = -6; xend = 0; xspan = [xstart xend];
cstart = 0; cstep = 0.1; cend = 45; C = cstart:cstep:cend;
intervals = 2 ^ 12;

U0 = zeros(1, length(C));
Ux0 = zeros(1, length(C));

for i = 1:length(C)
    init = asympt_left(params, C(i), xstart);
    [~, U] = f_solve(params, xspan, init, intervals);
    U0(i) = U(end, 1); Ux0(i) = U(end, 2);
end

figure('Position', [100 100 500 400]); hold on
xlabel('u(0)'); ylabel('u_x(0)')

% \gamma_- curve
plot([-U0(end:-1:1) U0], [-Ux0(end:-1:1) Ux0], '--', 'Color', 'black')

% \gamma_+ curve
plot([-U0(end:-1:1) U0], [Ux0(end:-1:1) -Ux0], 'Color', 'black')

end


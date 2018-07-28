%% Test (u(0), u_x(0)) - diagrams plotting
% Compare with Kolesnikova bachelor work
omega = 0; alpha = 1; params = [omega alpha];
diagram(params)

omega = 0; alpha = 8; params = [omega alpha];
diagram(params)

% Compare with G. L. Alfimov and D. A. Zezyulin
%   Nonlinear modes for the Gross-Pitaevskii equation -- a demonstrative computation approach

oemga = 0; alpha = 0; params = [omega, alpha];
diagram(params)
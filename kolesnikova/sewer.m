function [ X, U ] = sewer( params, cleft_init, cright_init )
% TODO: add description
% TODO: add internal parameters (like `intervals`, `tolx` and so on)
%
% INPUT:
%

% INTERNAL PARAMETERS:
tolx = 1e-9;
options = optimset('TolX', tolx);

xleft = -8; xright = -xleft; xend = 0;
xspan_left = [xleft xend]; xspan_right = [xright xend];
intervals = 2 ^ 12;


% Prepare internal fucntion for `fminsearch` usage
shots_distance_params = @(C) shots_distance(params, xleft, C(1), C(2));
[cmin, ~] = fminsearch(shots_distance_params, [cleft_init, cright_init], options);

init_left = asympt_left(params, cmin(1), xleft);
init_right = asympt_right(params, cmin(2), xright);

[X_left, U_left] = f_solve(params, xspan_left, init_left, intervals);
[X_right, U_right] = f_solve(params, xspan_right, init_right, intervals);

X = [X_left; X_right(end:-1:1)];
U = [U_left; U_right(end:-1:1, :)];

end


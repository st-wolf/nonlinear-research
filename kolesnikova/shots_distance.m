function [ distance ] = shots_distance( params, xleft, cleft, cright )
% TODO
%
% INPUT:
%   

% INTERNAL PARAMETERS:
xright = -xleft; xend = 0;
xspan_left = [xleft xend]; xspan_right = [xright xend];
intervals = 2 ^ 12;

init_left = asympt_left(params, cleft, xleft);
init_right = asympt_right(params, cright, xright);

[~, U_left] = f_solve(params, xspan_left, init_left, intervals);
[~, U_right] = f_solve(params, xspan_right, init_right, intervals);
 
distance = norm(U_left(end, :) - U_right(end, :));

end


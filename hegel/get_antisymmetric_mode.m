function [ X, U ] = get_antisymmetric_mode( params, C, xspan )
% Calculate a symmetric mode using the parameter C of the asymptotic 'asympt'
%
% INPUT:
%

xstart = xspan(1);

init = asympt_left(params, C, xstart);

[X, U] = f_solve(params, xspan, init, 2 ^ 14);

X = [X; -X(end-1:-1:1)];
U = [
	U(:, 1), U(:, 2);
	-U(end-1:-1:1, 1), U(end-1:-1:1, 2)
];

end
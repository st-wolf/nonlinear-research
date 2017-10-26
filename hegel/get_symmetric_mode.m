function [ X, U ] = get_symmetric_mode( mex_solver_name, params, C, xspan )
% Calculate a symmetric mode using the parameter C of the asymptotic 'asympt'
%
% INPUT:
%

mex_solver = str2func(mex_solver_name);

xstart = xspan(1);

init = asympt_left(params, C, xstart);

[X, U] = mex_solver(params, xspan, init, 2 ^ 12);

X = [X; -X(end-1:-1:1)];
U = [
	U(:, 1), U(:, 2);
	U(end-1:-1:1, 1), -U(end-1:-1:1, 2)
];

end
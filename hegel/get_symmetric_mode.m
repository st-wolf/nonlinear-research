function [ X, U ] = get_symmetric_mode( params, C, xspan )
% Calculate a symmetric mode using the parameter C of the asymptotic 'asympt'
%
% INPUT:
%

xstart = xspan(1);

init = [ asympt(params, C, xstart), asympt_diff(params, C, xstart) ];

[X, U] = f_solve(params, xspan, init, 1024);

X = [X -X(end-1:-1:1)];
U = [
	U(:, 1), U(:, 2);
	U(end-1:-1:1, 1), -U(end-1:-1:1, 2)
];

end
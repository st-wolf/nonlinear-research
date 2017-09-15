function [ X, U ] = get_mode( params, C, xspan )

xleft = xspan(1);
xright = xspan(2);

intervals = 2^10;
eps = 1e-6;

init_left = [asympt_left(params, C, xleft)
	(asympt_left(params, C, xleft + eps) - asympt_left(params, C, xleft - eps)) / (2 * eps)];

[X_left, U_left] = f_solve(params, [xleft 0], init_left, intervals);

init_right = [asympt_right(params, C, xright)
	(asympt_right(params, C, xright + eps) - asympt_right(params, C, xright - eps)) / (2 * eps)];

[X_right, U_right] = f_solve(params, [xright 0], init_right, intervals);

X = [X_left; X_right(end-1:-1:1)];

U = [
	U_left(:, 1), U_left(:, 2);
	U_right(end-1:-1:1, 1), U_right(end-1:-1:1, 2)
];

end


function difference = asympt_diff( params, C, xspan )

xleft = xspan(1);
xright = xspan(2);

intervals = 2^10;
eps = 1e-6;

init_left = [asympt_left(params, C, xleft)
	(asympt_left(params, C, xleft + eps) - asympt_left(params, C, xleft - eps)) / (2 * eps)];

[~, U] = f_solve(params, [xleft 0], init_left, intervals);

u_left = U(end, 1);
ux_left = U(end, 2);

init_right = [asympt_right(params, C, xright)
	(asympt_right(params, C, xright + eps) - asympt_right(params, C, xright - eps)) / (2 * eps)];

[~, U] = f_solve(params, [xright 0], init_right, intervals);

u_right = U(end, 1);
ux_right = U(end, 2);

difference = norm([u_left ux_left] - [u_right ux_right]);

end


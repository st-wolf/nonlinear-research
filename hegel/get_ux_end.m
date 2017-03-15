function ux_end = get_ux_end( params, C, xspan )

xstart = xspan(1);
init = [asympt(params, C, xstart), asympt_diff(params, C, xstart)];

[~, U] = f_solve(params, xspan, init, 1024);

ux_end = U(end, 2);

end
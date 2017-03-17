function ux_end = get_ux_end( params, C, xspan )

xstart = xspan(1);
init = asympt_left(params, C, xstart);

[~, U] = f_solve(params, xspan, init, 2 ^ 12);

ux_end = U(end, 2);

end
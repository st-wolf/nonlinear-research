function u_end = get_u_end( params, C, xspan )

eps = 1e-6;

xstart = xspan(1);
init = asympt_left(params, C, xstart);

[~, U] = f_solve(params, xspan, init, 2 ^ 12);

u_end = U(end, 1);

end


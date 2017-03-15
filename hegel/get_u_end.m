function u_end = get_u_end( params, C, xspan )

xstart = xspan(1);
init = [asympt(params, C, xstart), asympt_diff(params, C, xstart)];

[~, U] = f_solve(params, xspan, init, 1024);

u_end = U(end, 1);

end


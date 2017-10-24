function ux_end = get_ux_end( mex_solver, params, C, xspan )

xstart = xspan(1);
init = asympt_left(params, C, xstart);

[~, U] = mex_solver(params, xspan, init, 2 ^ 14);

if isempty(U)
    ux_end = NaN;
else
    ux_end = U(end, 2);
end

end
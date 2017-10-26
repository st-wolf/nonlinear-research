function u_end = get_u_end( mex_solver_name, params, C, xspan )

xstart = xspan(1);
init = asympt_left(params, C, xstart);

mex_solver = str2func(mex_solver_name);
[~, U] = mex_solver(params, xspan, init, 2 ^ 14);

if isempty(U)
    u_end = NaN;
else
    u_end = U(end, 1);
end

end


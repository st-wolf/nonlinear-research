function [ phi_norm ] = compute_asymmetric_mode_norm( params )
%
%

xspan = [-pi 0];

c_asymmetric = get_asymmetric_mode_parameter(params, xspan);
[X, Phi] = get_asymmetric_mode(params, c_asymmetric, xspan);
phi_norm = get_solution_norm(X, Phi);

end



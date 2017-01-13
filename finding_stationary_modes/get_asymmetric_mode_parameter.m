function [ c_asymmetric ] = get_asymmetric_mode_parameter( params, xspan )
% 
%
% INPUT:
%	

get_phi_end_params = @(C) get_phi_end(params, C, xspan);

c0 = 0.01; cstep = 0.1;
phi_end_0 = get_phi_end_params(c0);

c1 = c0 + cstep;
phi_end_1 = get_phi_end_params(c1);

while sign(phi_end_0) == sign(phi_end_1)
	c0 = c1; phi_end_0 = phi_end_1;
	
	c1 = c0 + cstep;
	phi_end_1 = get_phi_end_params(c1);
end

eps = 1e-5;
c_asymmetric = dichotomy(get_phi_end_params, c0, c1, eps);

end
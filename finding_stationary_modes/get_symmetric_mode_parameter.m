function [ c_symmetric ] = get_symmetric_mode_parameter( params, xspan )
% 
%
% INPUT:
%	

get_phix_end_params = @(C) get_phix_end(params, C, xspan);

c0 = 0.01; cstep = 0.1;
phix_end_0 = get_phix_end_params(c0);

c1 = c0 + cstep;
phix_end_1 = get_phix_end_params(c1);

while sign(phix_end_0) == sign(phix_end_1)
	c0 = c1; phix_end_0 = phix_end_1;
	
	c1 = c0 + cstep;
	phix_end_1 = get_phix_end_params(c1);
end

eps = 1e-5;
c_symmetric = dichotomy(get_phix_end_params, c0, c1, eps);

end


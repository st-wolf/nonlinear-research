function [ X, Phi ] = get_mode_from_cspan( params, xspan, cspan )
%
%

get_phix_end_params = @(C) get_phix_end(params, C, xspan);

eps = 1e-5;
c = dichotomy(get_phix_end_params, cspan(1), cspan(2), eps);

[ X, Phi ] = get_symmetric_mode(params, c, xspan);

end


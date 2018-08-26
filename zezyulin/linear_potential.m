function potential = linear_potential( params, xgrid )
% TODO:
%
% INPUT:
%

[~, omega, ~] = parse_params(params);
potential = 0.5 * omega^2 * (xgrid .^ 2);

end


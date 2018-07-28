function potential = nonlinear_potential( params, xgrid )
% TODO:
%
% INPUT:
%

[~, alpha] = parse_params(params);
potential = -(1 + alpha * (xgrid .^ 2));

end


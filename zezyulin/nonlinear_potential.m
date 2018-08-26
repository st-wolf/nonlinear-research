function potential = nonlinear_potential( params, xgrid )
% TODO:
%
% INPUT:
%

[~, ~, A] = parse_params(params);
% potential = (1 + A * (tanh(xgrid) .^ 2));
potential = 1 + A * xgrid .^ 2;

end


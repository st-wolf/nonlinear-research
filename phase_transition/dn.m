function [ f ] = dn( x, lambda )
[~, ~, f] = ellipj(x, lambda);
end
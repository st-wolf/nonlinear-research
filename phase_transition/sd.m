function [ f ] = sd( x, lambda )
f = sn(x, lambda) ./ dn(x, lambda);
end


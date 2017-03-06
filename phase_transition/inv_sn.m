function [ f ] = inv_sn( x, lambda )
f = elliptic12(asin(x), lambda);
end


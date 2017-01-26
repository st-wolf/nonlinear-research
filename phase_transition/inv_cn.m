function [ f ] = inv_cn( x, lambda )
f = elliptic12(acos(x), lambda);
end


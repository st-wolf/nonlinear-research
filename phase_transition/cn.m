function [ f ] = cn( x, lambda )
[~, f, ~] = ellipj(x, lambda);
end
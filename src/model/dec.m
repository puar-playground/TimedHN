function [index] = dec(x)
% convert a binary vector to a decimal number index
n = length(x);

% the basis vector
base_vec = arrayfun(@(x) 2^(x-1), 1:1:n);
% the index (start from 1)
index = x * base_vec' + 1;
end

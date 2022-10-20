function [x] = dec_inv(index, n)
% convert a decimal number index to a binary vector
x = bitget(index-1, 1:n);
end
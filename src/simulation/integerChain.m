function [G] = integerChain(n)
% Generate adjacent matrix of the integer chain.
G = zeros(n);
G(1:end-1, 2:end) = eye(n-1);
end
function [f] = F_score(A, R_true)
% compute hamming distance between positive sparse matrix

p = precision(A, R_true);
r = recall(A, R_true);
f = 2 / (1/r + 1/p);

end
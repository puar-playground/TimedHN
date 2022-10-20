function [r] = recall(A, R_true)
% compute hamming distance between positive sparse matrix
A(A>0) = 1;
A(A<0) = 0;
R_true(R_true>0) = 1;
R_true(R_true<0) = 0;
A = A - diag(diag(A));
R_true = R_true - diag(diag(R_true));

correct = sum(sum(A .* R_true));
true_cnt = sum(sum(R_true));
r = correct / true_cnt;

end
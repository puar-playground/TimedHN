function [p] = precision(A, R_true)
% compute hamming distance between positive sparse matrix
A(A>0) = 1;
A(A<0) = 0;
R_true(R_true>0) = 1;
R_true(R_true<0) = 0;
A = A - diag(diag(A));
R_true = R_true - diag(diag(R_true));

correct = sum(sum(A .* R_true));
prediction_cnt = sum(sum(A));
p = correct / prediction_cnt;

end
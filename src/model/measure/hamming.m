function [d] = hamming(A, B)
% compute hamming distance between positive sparse matrix
A(A>0) = 1;
A(A<0) = 0;
B(B>0) = 1;
B(B<0) = 0;
A = A - diag(diag(A));
B = B - diag(diag(B));

diff = (A ~= B);
d = sum(sum(diff));

end
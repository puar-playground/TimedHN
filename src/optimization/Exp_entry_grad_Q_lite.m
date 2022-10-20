function [Q_grad] = Exp_entry_grad_Q_lite(t, Q, i, j)
% comuter the gradient of the i,j-th entry of the matrix exponential expm(t*Q)
% Which is wirtten as: d expm(t*Q)_{ij} / d Q

% directional matrix V, only (i,j) entry is 1, 0 elsewhere.
V = zeros(size(Q));
V(i, j) = 1;
n = size(Q, 1);
A = repmat(t*Q', 2, 2);
A(n+1:end, 1:n) = 0;
A(1:n, n+1:end) = V;
Q_grad_temp = t*expm(A);
Q_grad = Q_grad_temp(1:n, n+1:end);

end
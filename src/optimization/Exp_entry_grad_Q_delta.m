function [Q_grad] = Exp_entry_grad_Q_delta(t, Q, i, j)
% comuter the gradient of the i,j-th entry of the matrix exponential expm(t*Q)
% Which is wirtten as: d expm(t*Q)_{ij} / d Q

% directional matrix V, only (i,j) entry is 1, 0 elsewhere.
V = zeros(size(Q));
V(i, j) = 1;
n = size(Q, 1);
delta = 0.0001;
Q_grad = (expm(t*(Q'+delta*V)) - expm(t*Q')) / delta;

end
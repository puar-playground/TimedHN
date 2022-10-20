function [Q_grad] = Exp_entry_grad_Q(t, Q, i, j)
% comuter the gradient of the i,j-th entry of the matrix exponential expm(t*Q)
% Which is wirtten as: d expm(t*Q)_{ij} / d Q

% directional matrix V, only (i,j) entry is 1, 0 elsewhere.
V = zeros(size(Q));
V(i, j) = 1;
% the integral from 0 to t
d_Q_grad = @(tau) expm((t-tau)*Q')*V*expm(tau*Q');
% the gradient
Q_grad = integral(d_Q_grad, 0, t, 'ArrayValued', true);

% Q_binary = zeros(size(Q));
% Q_binary(Q~=0) = 1;
% 
% Q_grad = sparse(Q_binary.*Q_grad);

end
function [R_grad, likelihood] = sample_grad(observation, t, Q)
% compute the gradient w.r.t R, using an observation. 

% get number of events.
n_event = length(observation);

% Permutation and compute sub-P_t
[Q_sub, ~, T_col] = sample_permute_Q(observation, Q);
P_t_sub = expm(t*Q_sub);
likelihood = P_t_sub(1, end);

%  Exponential Gradient wrt Q
Q_grad = Exp_entry_grad_Q(t, Q_sub, 1, size(Q_sub, 2));
Q_grad(size(Q, 1), size(Q, 2)) = 0;
Q_grad_global = T_col * Q_grad * T_col';

% Q Gradient wrt R
R_grad = Q_grad_R(Q_grad_global, n_event) / likelihood;

end
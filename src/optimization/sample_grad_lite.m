function [R_grad, t_grad, likelihood] = sample_grad_lite(observation, t, R)
%
k = sum(observation);
Q_sub = Q_sub_lite(observation, R);
P_t_sub = expm(t*Q_sub);

% compute the dP/dt = Qexp(tQ)
t_grad_all = Q_sub * P_t_sub;
t_grad = t_grad_all(1, 2^k);

likelihood = P_t_sub(1, end);
if k > 7
    Q_grad_sub = Exp_entry_grad_Q_delta(t, Q_sub, 1, 2^k);
else
    Q_grad_sub = Exp_entry_grad_Q_lite(t, Q_sub, 1, 2^k);
end
[rows, cols] = Q_sparse_position(sum(observation));
% [rows, cols] = find(Q_sub);

R_grad = lite_grad_R(observation, Q_grad_sub, rows, cols);

end


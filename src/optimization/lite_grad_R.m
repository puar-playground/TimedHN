function [R_grad] = lite_grad_R(observation, Q_grad_sub, rows, cols)
% compute gradient w.r.t R using the matrix exponential in the subspace

% n is the number of genes.
n_event = length(observation);
% extract all the mutated genes.
mut_genes = find(observation == 1);
% count the number of mutated genes.
n_mut = sum(observation);

% it has 2^n_mut different mutation states.
n_state_sub = size(Q_grad_sub, 1);
% the absorbing states (fully mutated) is the last state.
final_state_sub = 2^n_mut - 1;
final_state = 2^n_event - 1;

% empty matrix to store dL / dR
R_grad = zeros(n_event);

for i=1:length(rows)
    % possible transition index pair
    r = rows(i);
    c = cols(i);
    
    % bit vectors in subspace
    rb_sub = dec_inv(r, n_mut);
    cb_sub = dec_inv(c, n_mut);
    
    % gradient of the entry
    mute_after = mut_genes(cb_sub==1);
    to_mute = mut_genes(cb_sub - rb_sub==1);
    entry_grad = zeros(n_event);
    entry_grad(mute_after, to_mute) = 1;

    % update R_grad
    R_grad = R_grad + Q_grad_sub(r, c) * entry_grad;
    
end

for i=1:n_state_sub
    bit_sub = dec_inv(i, n_mut);
    state_vec = zeros(1, n_event);
    state_vec(mut_genes(bit_sub==1)) = 1;
    diag_index = dec(state_vec);

    entry_grad = diag_Q_grad_R(diag_index, n_event);
    R_grad = R_grad + Q_grad_sub(i, i) * entry_grad;
    
end

end




function [entry_grad] = diag_Q_grad_R(diag_index, n_event)
% compute the accumulated Q_grad_R on diagonal entry.

final_state = 2^n_event - 1;
rb = dec_inv(diag_index, n_event);
entry_grad = zeros(n_event);
to_genes = find(bitget(bitxor(final_state, diag_index-1) , 1:n_event) == 1);
for i=1:length(to_genes)
    entry_grad_temp = zeros(n_event);
    cb = rb;
    cb(to_genes(i)) = 1;
    entry_grad_temp((cb == 1), (cb - rb == 1)) = 1;
    entry_grad = entry_grad + entry_grad_temp;
end

entry_grad = -entry_grad;

end
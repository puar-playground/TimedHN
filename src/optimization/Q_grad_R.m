function [R_grad] = Q_grad_R(Q_grad_global, n_event)
% compute dQ /dR

R_grad = zeros(n_event);
final_state = 2^n_event - 1;

[row_i,col_j,vs] = find(Q_grad_global);

for i=1:length(vs)

    row = row_i(i);
    col = col_j(i);
    rb = dec_inv(row, n_event);
    cb = dec_inv(col, n_event);

    % check transition validity
    if (max(cb - rb) ~= sum(cb - rb)) || min(cb - rb) ~= 0
        continue
    elseif (max(cb - rb) == sum(cb - rb)) && max(cb - rb) == 1
        entry_grad = zeros(n_event);
        entry_grad((cb == 1), (cb - rb == 1)) = 1;
        R_grad = R_grad + vs(i) * entry_grad;
    
    elseif (max(cb - rb) == sum(cb - rb)) && max(cb - rb) == 0
        % row == col
        entry_grad = diag_Q_grad_R(row, n_event);
        R_grad = R_grad + vs(i) * entry_grad;
    end
        

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
function [rows, cols] = Q_sparse_position(n_event)

% it has 2^n different mutation states.
n_state = 2^n_event;
% the absorbing states (fully mutated) is the last state.
final_state = 2^n_event - 1;

% Initialize row, column, value arrays for the sparse representation of Q.
% Q has n*2^(n-1) number of values.
rows = int32.empty(n_event*2^(n_event-1), 0);
cols = int32.empty(n_event*2^(n_event-1), 0);

% Iteratively fill the matrix.
cnt = 1;
% traverse all the states. (the decimal values of the binary values)
for i=0:n_state - 1
    % extract all the unmutated genes' indices.
    to_genes = find(bitget(bitxor(final_state, i) , 1:n_event) == 1);

    % traverse all the unmutated genes for the state i.
    for j=1:length(to_genes)
        % the i+1-th states
        rows(cnt) = i+1;
        % the index of the target state.
        cols(cnt) = bitset(i, to_genes(j))+1;
        % increment saving index.
        cnt = cnt + 1;
    end

end

end
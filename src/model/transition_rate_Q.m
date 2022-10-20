function [Q] = transition_rate_Q(R)

n = size(R, 1);
if n > 12
    fprintf('compute Q in parallel\n')
    Q = transition_rate_Q_parallel(R);
else
    Q = transition_rate_Q_sequential(R);
end

end

function [Q] = transition_rate_Q_sequential(R)
% Fill Q using R
% All mutation states are represented by bit arrays and sorted by their value
% as binary-numbers.

% n is the number of genes.
n = size(R, 1);
% it has 2^n different mutation states.
n_state = 2^n;
% the absorbing states (fully mutated) is the last state.
final_state = 2^n - 1;

% Initialize row, column, value arrays for the sparse representation of Q.
% Q has n*2^(n-1) number of values.
rows = int32.empty(n*2^(n-1), 0);
cols = int32.empty(n*2^(n-1), 0);
values = double.empty(n*2^(n-1), 0);

% Iteratively fill the matrix.
cnt = 1;
% traverse all the states. (the decimal values of the binary values)
for i=0:n_state - 1
    % extract all the unmutated genes' indices.
    to_genes = find(bitget(bitxor(final_state, i) , 1:n) == 1);
    % diag is a variable recording the i-th diagonal value of Q.
    diag_v = 0;

    % the bit array of the i-th state. 
    % (the first number is the first gene)
    row_state = bitget(i, 1:n);

    % traverse all the unmutated genes for the state i.
    for j=1:length(to_genes)
        % the i+1-th states
        rows(cnt) = i+1;
        % the index of the target state.
        cols(cnt) = bitset(i, to_genes(j))+1;
        

        % compute the collective transition rate. (competing exponentials)
        v = row_state * R(:, to_genes(j)) + R(to_genes(j), to_genes(j));
        values(cnt) = v;

        % update diagonal values.
        diag_v = diag_v + v;

        % increment saving index.
        cnt = cnt + 1;
    end

    % save diagonal to row, column, value arrays of the sparse Q.
    rows(cnt) = i+1;
    cols(cnt) = i+1;
    values(cnt) = -diag_v;
    % increment saving index
    cnt = cnt + 1;
end

% Convert row, column, value arrays to sparse matrix.
Q = sparse(rows, cols, values, 2^n, 2^n);

end


function [Q] = transition_rate_Q_parallel(R)
% Fill Q using R
% All mutation states are represented by bit arrays and sorted by their value
% as binary-numbers.

% n is the number of genes.
n = size(R, 1);
n_sub_state = 2^(n-1);

% initiallize sparse matrix value arrays
rows = zeros(n, 2^(n-1));
cols = zeros(n, 2^(n-1));
values = zeros(n, 2^(n-1));

parfor mut=1:n
    % add mut-th mutation
    base_vec = arrayfun(@(x) 2^(x-1), [1:mut-1, mut+1:n]);
    rate_array = R(:, mut);

    rows_temp = int32.empty(2^(n-1), 0);
    cols_temp = int32.empty(2^(n-1), 0);
    values_temp = double.empty(2^(n-1), 0);

    for remain_dec=1:n_sub_state
        % traverse the 2^(n-1) possible states for the remaining (n-1) bits 
        remain_vec = bitget(remain_dec-1, 1:n-1);
        
        % get the row index
        row_dec = remain_vec * base_vec' + 1;
        % get the column index
        col_dec = row_dec + 2^(mut - 1);
        % get the column state vector
        col_vec = bitget(col_dec-1, 1:n);
        % compute Q(row, col)
        v = col_vec * rate_array;

        % update the sparse arrays
        rows_temp(remain_dec) = row_dec;
        cols_temp(remain_dec) = col_dec;
        values_temp(remain_dec) = v;
    end

    rows(mut, :) = rows_temp;
    cols(mut, :) = cols_temp;
    values(mut, :) = values_temp;

end

% Convert row, column, value arrays to sparse matrix.
Q = sparse(rows(:), cols(:), values(:), 2^n, 2^n);
% Q(2^n, 2^n) = 0;

% updates the diagonal
Q = Q - diag(sum(Q, 2));

end
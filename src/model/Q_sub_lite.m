function [Q_sub] = Q_sub_lite(observation, R)
% this is a lite function to compute Q_sub from a given observation 
% without computing the full Q matrix.

% n is the number of genes.
n_event = length(observation);
% extract all the mutated genes.
mut_genes = find(observation == 1);
% count the number of mutated genes.
n_mut = sum(observation);

% it has 2^n_mut different mutation states.
n_state_sub = 2^n_mut;
% the absorbing states (fully mutated) is the last state.
final_state_sub = 2^n_mut - 1;
final_state = 2^n_event - 1;


% spontaneous
R_spontaneous = diag(R);


% Initialize row, column, value arrays for the sparse representation of Q.
% Q has n*2^(n-1) number of values.
rows = int32.empty(n_mut*2^(n_mut-1), 0);
cols = int32.empty(n_mut*2^(n_mut-1), 0);
values = double.empty(n_mut*2^(n_mut-1), 0);

% Iteratively fill the matrix.
cnt = 1;
% traverse all the states. (the decimal values of the binary values)
for i=0:n_state_sub - 1
    % extract indices of all the unmutated genes in subspace.
    to_genes_sub = find(bitget(bitxor(final_state_sub, i) , 1:n_mut) == 1);
    
    row_state = zeros(1, n_event);
    row_state(1, mut_genes(bitget(i, 1:n_mut) == 1)) = 1;

    % traverse all the unmutated genes for the state i.
    for j=1:length(to_genes_sub)
        % the i+1-th states
        rows(cnt) = i+1;
        % the index of the target state.
        cols(cnt) = bitset(i, to_genes_sub(j))+1;
        % the bit array of the i-th state. 
        % (the first number is the first gene)
        

        mut_gene = mut_genes(to_genes_sub(j));
        

        % compute the collective transition rate. (competing exponentials)
        v = row_state * R(:, mut_gene) + R(mut_gene, mut_gene);
        values(cnt) = v;

        % increment saving index.
        cnt = cnt + 1;
    end

    % extract indices of all the unmutated genes globally.
    to_genes = find(1 - row_state == 1);
    % diagnonal value
    rows(cnt) = i+1;
    cols(cnt) = i+1;
    values(cnt) = - sum(row_state * R(:, to_genes)) - sum(R_spontaneous(to_genes));
    cnt = cnt + 1;

end

% Convert row, column, value arrays to sparse matrix.
Q_sub = sparse(rows, cols, values);

end
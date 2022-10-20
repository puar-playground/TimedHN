function [Q_sub, state_index, T_col] = sample_permute_Q(observation, Q)

% This is a function compute the likelihood of a sample.
% Q_sub is the transition probability matrix in the sub-space

% n is the number of genes.
n = length(observation);
% It has 2^n different mutation states.
n_state = 2^n;
% Initialize an array to record the swaps.
state_index = 1:n_state;
% % the index of the observed state.
% absorb_index = bin2dec(sprintf('%d',observation(n:-1:1))) + 1;

% extract all the mutated genes.
mut_genes = find(observation == 1);
% count the number of mutated genes.
n_mut = sum(observation);
% a exponential bases array map indices of Q to the indices of the subspace of Q.
% the subspace is supported only by the mutated genes.
base_vec = arrayfun(@(x) 2^(x-1), mut_genes);

% The subspace (the first 2^(n_mut) diagonal block).
% Traverse all the sub-states.
for sub_index=1:2^(n_mut) - 1
    % Compute the index on the Q matrix.
    glob_index = bitget(sub_index, 1:n_mut) * base_vec';
    % Record the swap using the array "state_index"
    state_index(glob_index+1) = sub_index+1;
    state_index(sub_index+1) = glob_index+1;
%     % Print all the swaps if needed
%     fprintf('state %i swap to state %i\n', old_index+1, i+1)
end

% Initialize an identity matrix.
T = speye(n_state);
% Make it a permutation matrix using the swaps
T_col = sparse(T(:, state_index));
% Modify the order of states
Q = sparse(T_col' * Q * T_col);
% Get the subspace (the first 2^(n_mut) diagonal block)
% the transition from zero state to the observed state is stored in the
% first row and the last column of Q_sub(1, 2^(n_mut))
Q_sub = Q(1:2^(n_mut), 1:2^(n_mut));

state_index = state_index(1:2^(n_mut));

end
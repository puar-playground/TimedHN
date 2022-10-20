function [adj_matrix] = randDAG(n_event, n_edge)
% rand DAG 
adj_matrix = zeros(n_event);
n_use = n_event - 1;
edge_pool = randsample(1:n_use*(n_use+1)/2, n_edge);
for i = 1:n_edge
    edge_index = edge_pool(i);
    c = ceil((sqrt(1+ 4 * 2 * edge_index) - 1) / 2);
    r = edge_index - ((c-1)^2 + c-1) / 2;
    adj_matrix(r, c+1) = 1;
end
spontaneous = zeros(1, n_event);
spontaneous(1, sum(adj_matrix)==0) = 1;
adj_matrix = adj_matrix + diag(spontaneous);

end
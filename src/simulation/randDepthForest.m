function [adj_matrix] = randDepthForest(n_event)
% rand DAG with log2(n) depth with max number of trees
n_layer = round(log2(n_event)) + 1;
max_tree = n_layer - 1;

if max_tree>2
    n_tree = randsample(2:max_tree, 1);
else
    n_tree = 2;
end

layers = cell(1,n_layer);
layers{1} = 1:n_tree;

depth_splite = sort(randsample(n_tree+1:n_event-2, n_layer-2));
depth_splite = [n_tree, depth_splite, n_event];

for l=2:n_layer
    layers{l} = depth_splite(l-1)+1:depth_splite(l);
end

adj_matrix = zeros(n_event);

for d=n_layer:-1:2
    child_layer_node = layers{d};
    parent_layer_node = layers{d-1};

    for j=1:length(child_layer_node)
        child = child_layer_node(j);
        
        if length(parent_layer_node)>1
            parent = randsample(parent_layer_node, 1);
        else
            parent = parent_layer_node(1);
        end
        adj_matrix(parent, child) = 1;
    end

end

root_layer = layers{1};
for i=1:length(root_layer)
    root = root_layer(i);
    adj_matrix(root, root) = 1;
end

end

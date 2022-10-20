function [] = show_result_single(R_best, gene_name)
 
% Plot the inferred result and the truth.
figure()
% R_best = R_best - diag(diag(R_best));
dependent = sum(R_best) + sum(R_best');
dependent_index = dependent > 0;
R_best = R_best(dependent_index, dependent_index);
gene_name = gene_name(dependent_index);

g = digraph(R_best);
g.Nodes.Name = gene_name';
LWidths = 2*g.Edges.Weight/max(g.Edges.Weight);
plot(g,'Layout','layered', 'EdgeLabel',g.Edges.Weight,'LineWidth',LWidths, ...
    'Marker', 'o', 'NodeColor','k', 'MarkerSize', 10, 'LineStyle', '-', ...
    'ArrowSize', 12, 'EdgeColor', 'b')

end
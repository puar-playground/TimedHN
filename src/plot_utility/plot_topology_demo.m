clc
clear
close all
addpath(genpath(pwd))
rng(1)

%%
n_event = 15;
structure = 'Forest_new';
run = 2;

figure('Position', [10 10 500 540])
subplot(1, 1, 1, 'Position', [0.02 0.12 0.96 0.8])
R_true = readtable(strcat(pwd, '/simulation_data/', structure, '/R_', num2str(run), '.csv'));
R_true = R_true{:, :};
sp = diag(R_true);
R_true = R_true - diag(diag(R_true))
% R_true = R_true - diag(0.1*ones(1, n_event));

g = digraph(R_true);
ms = (sp + 0.1 * ones(15, 1))*10 + 10;
h = plot(g,'Layout','layered','LineWidth',4*g.Edges.Weight, ...
    'Marker', 'o', 'NodeColor',[0.4,0.4,0.4], 'MarkerSize', ms, 'LineStyle', '-', ...
    'ArrowSize', 15, 'EdgeColor', 'k', 'ArrowPosition', 0.9)

labelnode(h, 1:n_event, '')
box on;
ax = gca;
ax.LineWidth = 4;
title(strcat('Branching'), 'FontSize',35)
saveas(gcf, strcat(pwd, '/plot/', num2str(n_event), '_', structure, '_topology.png'))
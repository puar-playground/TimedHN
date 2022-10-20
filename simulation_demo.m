clc
clear
close all
addpath(genpath(pwd))
% rng(1)

%%
n_event = 15;
max_mut = 10;
n_sample = 1000;
structure = 'DAG';

R = randDepthForest(n_event); % for forest structure
% R = randDAG(n_event, floor(1.5*n_event)); % for DAG structure
R = R + diag(0.01*ones(1, n_event));

% save generated graph adjacency matrix R
R_table = array2table(R);
header = string(1:n_event);
header = strcat('gene_', header);
R_table.Properties.VariableNames = header;
writetable(R_table, strcat(pwd, '/data/simiulation_R.csv'))

% Sampling synthetic data
samples = CTMCrnd(R, n_event, n_sample, max_mut, 'uniform');

% save synthetic data as csv file
snapshort = array2table(samples);
snapshort.Properties.VariableNames(1) = {'Time'};
header = string(1:n_event);
header = strcat('gene_', header);
snapshort.Properties.VariableNames(2:end) = header;
writetable(snapshort, strcat(pwd, '/data/simulation_data.csv'))

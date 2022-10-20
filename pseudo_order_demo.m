clc
clear
close all
addpath(genpath(pwd))
rng(0)

%%
R = readtable(strcat(pwd, '/result/hazard_network_R.csv'));
R = R{:, :};
snapshot = readtable(strcat(pwd, '/data/luminal_breast_cancer.csv'));
patients = string(snapshot{:, 1});
genes = string(snapshot.Properties.VariableNames);
genes = genes(2:end);
data = snapshot{:, 2:end};
states = zeros(1, size(data, 1));
for r=1:size(data, 1)
    states(1, r) = dec(data(r,:));
end

data_unique = unique(data, 'rows');
n_sample = size(data_unique, 1);

%%
order_strs = [];
patient_lists = cell(1, n_sample);
t_list = [];
table_height = 0;
tic
for i=1:n_sample
    observe = data_unique(i, :);
    state = dec(observe);
    patient_sub = patients(states==state, 1)
    table_height = max(table_height, length(patient_sub));
    t_list = [t_list, TimeExpectation(observe, R)];
    [max_like_order, max_likelihood] = max_likelihood_order(observe, R);
    order_strs = [order_strs, 'R --> ' + strjoin(genes(max_like_order), ' --> ')];
    patient_lists{1, i} = patient_sub;
end
toc

%% sort profile by t
data_unique_t = [t_list', data_unique];
[~, sorted_order] = sort(t_list);
data_unique_t = data_unique_t(sorted_order, :);
t_sorted_profile = array2table(data_unique_t);
t_sorted_profile.Properties.VariableNames = ['time_expectation p(t|x)', genes]
writetable(t_sorted_profile, strcat(pwd, '/result/pseudo_time_patients.csv'))

%% write orders and patients to csv
M = containers.Map(order_strs,patient_lists);
sz = [table_height+1, n_sample];
v_type = cell(1, n_sample);
v_type(:) = {"string"};
v_type = string(v_type);
order_table = table('Size',sz, 'VariableTypes', v_type, 'VariableNames', order_strs);
for i = 1:n_sample
    order_table(2:length(patient_lists{i}) + 1, i) = cellstr(patient_lists{i});
    order_table(1, i) = {num2str(t_list(i))};
end

order_strs_sorted = order_strs(1, sorted_order);
order_table = order_table(:, order_strs_sorted);

writetable(order_table, strcat(pwd, '/result/pseudo_time_events.csv'))
    






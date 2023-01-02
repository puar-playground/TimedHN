clc
clear
close all
addpath(genpath(pwd))
rng(0)

%% load data
snapshot = readtable(strcat(pwd, '/data/luminal_breast_cancer.csv'));
% load data profiles
data = snapshot{:, 2:end};
% show unique profiles
plot_binary(unique(data,'rows'));
gene_name = snapshot.Properties.VariableNames(2:end);
n_sample = size(snapshot, 1);
n_event = size(snapshot, 2) - 1;
patience_max = 1;
lambda = 1e-2;
lr_0 = 1e-4; % initial learning rate
decay_rate = 1e-4;
%% Initialize CTMC model
% use patrick suppes' conditions of probabilistic causation to initialize R
% (prima facie causality used in CAPRI)
R = 1 * suppes_init(data) + diag(0.1 * ones(1, n_event));
% scale sum of R down to 5, for faster convergence.
R = R * 5 / sum(sum(R));

% initialize times using conditional expectation E(t|x), such that the
% times and R are compatible. (optimization more stable)
times = zeros(1, n_sample);
parfor row=1:n_sample
    observation = data(row, :);
    times(row) = TimeExpectation(observation, R);
end
t_budget = sum(times);
R_best = R;

%% optimize
% initialize objective function
max_obj = -inf;
obj = -inf;
% patience for convergence
patience = 0;

tStart = tic;
% gradient update for 10000 iterations
while true

    lr = (1 / (1 + decay_rate * iter)) * lr_0;
    [R_grad_total, t_grad_total, log_likelihood] = batch_grad(data, times, R);

    % compute objective
    L1_global = sum(sum(R));
    obj = log_likelihood - lambda * L1_global;
    
    % update R_best
    if obj - max_obj > 1*10^(-6)
        max_obj = obj;
        patience = 0;
        R_best = R;
    elseif obj - max_obj > 0
        patience = patience + 1;
    end
    
    % print iter index
    fprintf(['Iteration %i, obj: %f, last log-like: %f, L1_global: %f' ...
        ' cost %f second\n'], iter, obj, log_likelihood, L1_global, toc(tStart))
    
    % converged or not
    if patience >= patience_max
        break
    end

    % L1 gradient
    L1_grad = ones(n_event);
    
    % gradient update R and times
    grad_total = R_grad_total - lambda * L1_grad;
    R = R + lr * grad_total;
    % keep R positive
    R = max(0, R);
    % gradient update t
    times = times + lr * t_grad_total;
    % map time to positive
    times = max(0, times); 
    % map time to the constraint space
    t_sum = sum(times);
    times = times * (t_budget / t_sum);
    
    iter = iter + 1;
end

%% Plot the result after thresholding
spontaneous_r = diag(diag(R_best));
intergraph = R_best - spontaneous_r;
intergraph(intergraph <= 0.01) = 0;
R_best_th = intergraph + spontaneous_r;
% plot R after thresholding
show_result_single(R_best_th, gene_name);

%% save R result
R_table = array2table(R_best_th);
R_table.Properties.VariableNames = gene_name;
writetable(R_table, strcat(pwd, '/result/hazard_network_R.csv'))





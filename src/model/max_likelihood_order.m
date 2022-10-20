function [max_like_order, max_likelihood] = max_likelihood_order(observation, R)
% This function find the order of genes in a mutation profile that
% accumulates with the maximum likelihhod.

% number of all genes
n_gene = length(observation);
% the spontaneous rates
spontaneous = diag(R)';
% find which gene is going to mutate in the observation
mut = find(observation==1);
n_mut = sum(observation);
all_perms = perms(mut);
% initialize the max likelihood and order
max_like_order = mut;
max_likelihood = 0;
for i=1:size(all_perms, 1)
    likelihood = 1;
    order = all_perms(i, :);
    % mut_temp is a list of already-mutated genes
    mut_temp = [];
    for j=1:n_mut
        % compute rates to all un-mutated genes
        to_mut = setdiff(1:n_gene, mut_temp);
        R_sub = R(mut_temp, to_mut);
        stimulate = sum(R_sub);
        rates_to_all = stimulate + spontaneous(to_mut);
        
        % find the index of the next mut gene in the to_mut list
        new_mut = order(j);
        index_new = (to_mut == new_mut);
        p = rates_to_all(index_new) / sum(rates_to_all);
        likelihood = likelihood * p;

        % update the mut_temp list
        mut_temp = [mut_temp, new_mut]; 
    end
    if likelihood > max_likelihood
        max_likelihood = likelihood;
        max_like_order = order;
    end

end

end
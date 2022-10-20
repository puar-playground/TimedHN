function [R_grad_total, t_grad_total, log_likelihood] = batch_grad(data, tumor_age, R)
% compute average gradient and log-likelihood use a set of datas

n_sample = size(data, 1);
n_event = size(data, 2);
log_likelihood_array = zeros(1, n_sample);
R_grad_total_array = zeros(n_event, n_event, n_sample);
t_grad_total = zeros(1, n_sample);


% all positive samples
parfor s=1:n_sample
    observation = data(s, :);
    % gradient w.r.t R on a single observation
    [R_grad, t_grad, likelihood] = sample_grad_lite(observation, tumor_age(s), R);
    t_grad_total(1, s) = t_grad / (n_sample * likelihood);
    R_grad = R_grad / likelihood;
    % accumulate log_likelihood and gradient
    log_likelihood_array(1, s) = log(likelihood);
    R_grad_total_array(:, :, s) = R_grad;
end

% average gradient
R_grad_total = mean(R_grad_total_array, 3);
log_likelihood = mean(log_likelihood_array, 2);

end
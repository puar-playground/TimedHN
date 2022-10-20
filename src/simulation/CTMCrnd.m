function [samples] = CTMCrnd(R, n_event, n_sample, max_mut, time_mod)
% batch sampling from a CTMC gouverned by Q matrix.

% Compute Q
tStart = tic;
Q = transition_rate_Q(R);
fprintf('Transition rate Q filling cost %f second\n', toc(tStart))

% initialize samples array for storage.
samples = zeros(n_sample, n_event + 1);

parfor s=1:n_sample
    % draw sample mulitiple times
    [tumor_state, tumor_age] = CTMCrnd_single(Q, n_event, max_mut, time_mod);
    % store in array.
    samples(s, :) = [tumor_age, tumor_state];
end

samples = sortrows(samples, 1);

end


function [tumor_state, tumor_age] = CTMCrnd_single(Q, n_event, max_mut, time_mod)
% point sampling from a CTMC gouverned by Q matrix.

t_total = 0;
transition_cnt = 0;
state_now = 1;
max_transition_cnt = min(max_mut, n_event);
t_history = double.empty(max_transition_cnt + 1, 0);
t_history(1) = 0;
state_history = double.empty(max_transition_cnt + 1, 0);
state_history(1) = 1;

% simulate transition to max step.
while transition_cnt < max_transition_cnt
    % the bit vector of the current state.
    state_now_bit = dec_inv(state_now, 3);
    % total transition rate.
    r = -Q(state_now, state_now);
    % check if it is an absorbing state. 
    if r == 0
        break
    end

    % wait t before next state transition.
    t = exprnd(1/r);
    % transition distribution (probability of transition to each reachable state)
    [~, cols, vs] = find(Q(state_now, state_now+1:end));
    transition_distribution_sparse = vs / r;
    
    % transition_distribution = full(Q(state_now, state_now+1:end)) / r;
    
    % sampling the next state.
    state_next = state_now + cols((mnrnd(1, transition_distribution_sparse)==1));
    % state_next = state_now + find(mnrnd(1, transition_distribution)==1);

    % update transition count.
    transition_cnt = transition_cnt + 1;

    % update t and record history
    t_total = t_total + t;
    t_history(transition_cnt + 1) = t_total;
    state_history(transition_cnt + 1) = state_next;
    
    % update current state
    state_now = state_next;
end

if strcmp(time_mod,'uniform')
    % draw a uniformly distributed tumor age.
    tumor_age = 1 * t_total * rand;
elseif strcmp(time_mod,'exp')
    % draw a exponentially distributed tumor age
    tumor_age = min(exprnd(2), t_total);
end

sampled_index = find(t_history < tumor_age, 1, 'last');
tumor_state = dec_inv(state_history(sampled_index), n_event);

end
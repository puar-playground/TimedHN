function [R_init] = suppes_init(data)
% Initialize R matrix using the weight of the CAPRESE method
n_sample = size(data, 1);
n_event = size(data, 2);

%% precompute some probabilities

% the joint probability P(i, j)
P_joint = data' * data / n_sample;
% the probability of each event P(i)
P = diag(P_joint);
% the conditional probability of P(j|i)
P_i_j = P_joint ./ repmat(P, 1, n_event);

% the joint probability P(no i, j)
P_joint_ni = (1 - data)' * data / n_sample;
% the probability of each event P(no i) = 1 - P(i)
nP = 1 - diag(P_joint);
% the conditional probability of P(j| no i)
P_ni_j = P_joint_ni ./ repmat(nP, 1, n_event);


% my filter
contained = zeros(n_event);
contained(P_i_j==1) = 1;
contain = zeros(n_event);
contain(P_i_j'~=1) = 1;
my_filter = 1 - contain .* contained;



%% compute CAPRESE weight
R_init = zeros(n_event);
for i=1:n_event
    for j=1:n_event
        if i==j
            continue
        else
            R_init(i, j) = (P_joint(i, i) > P_joint(j, j));
        end
    end
end
R_init = R_init .* (P_i_j > P_ni_j);
%% initialize diagonal spontaneous
one_step_states = data(sum(data, 2)==1, :);
spontaneous = sign(abs(sum(one_step_states, 1)));
R_init = R_init + diag(spontaneous);
R_init = max(R_init, 0);

end
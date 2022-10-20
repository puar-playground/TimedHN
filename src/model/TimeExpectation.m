function [et] = TimeExpectation(observation, R)
% compute the Mathematical expectation of the time of obersving a state

% check observation fullness.
if sum(observation) == length(observation)
    % the last state is absorbing state. The time expectation is infinity.
    et = 10^8;
else
    [Q_sub] = Q_sub_lite(observation, R);
    % int_P_sub = @(tau) tau * expm(tau*Q_sub);
    % Et_sub = integral(int_P_sub, 0, 100, 'ArrayValued', true)
    
    % compute the integral of: t*exp(Qt), for t from 0 to +infinty.
    Et_sub = full(inv(Q_sub * Q_sub));
    % normalize the probability by dividing the integral of: exp(Qt), for t from 0 to +infinty.
    Et_normalized_sub = Et_sub ./ (-inv(Q_sub));
    % time expectation of the observation is in the first row last column.
    et = Et_normalized_sub(1, end);
end

end
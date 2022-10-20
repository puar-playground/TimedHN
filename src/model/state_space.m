function [state_space] = state_space(n_event)
%UNTITLED2 Summary of this function goes here
state_space = zeros(1 + n_event, 2^n_event);
for s=1:2^n_event
    state_space(2:end, s) = dec_inv(s, n_event);
end

end
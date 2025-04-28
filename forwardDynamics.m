function twist_dot = forwardDynamics(twist, t, w, I)
%FORWARDDYNAMICS Summary of this function goes here
%   Detailed explanation goes here
twist_dot = I \ (w(t)-ad(twist)'*I*twist);
end


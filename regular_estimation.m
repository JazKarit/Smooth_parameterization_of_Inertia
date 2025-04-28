function [pi_, iterations, res_norm] = regular_estimation(acc,twist_avg,w,t)
%REGULAR_ESTIMATION Summary of this function goes here
%   Detailed explanation goes here
    pi_init = [1 0 0 0 1 1 1 0 0 0]';
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxIterations',50);
    [pi_, res_norm, b, c, output]  = lsqnonlin(@(pi_) getResidual(acc,twist_avg,w,t,pi_),pi_init,[],[],options);
    iterations = output.iterations;
end


function [pi_, iterations, res_norm] = regular_estimation(acc,twist_avg,w,t)
%REGULAR_ESTIMATION Summary of this function goes here
%   Detailed explanation goes here
    function [c,ceq] = mycon(x)
        c = 0*x;     % Compute nonlinear inequalities at x.
        ceq = 0*x;   % Compute nonlinear equalities at x
    end

    pi_init = [1 0 0 0 1 1 1 0 0 0]';
    A = [-1 0 0 0 0 0 0 0 0 0;];
    b = [0];
    Aeq = [0 0 0 0 0 0 0 0 0 0;];
    Beq = [0];
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxIterations',50);
    %options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIterations',50);
    [pi_, res_norm, b, c, output]  = lsqnonlin(@(pi_) getResidual(acc,twist_avg,w,t,pi_),pi_init,[],[],A,b,Aeq,Beq,@mycon,options);
    iterations = output.iterations;
end


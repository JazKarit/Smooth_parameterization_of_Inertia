function [pi_,iterations,res_norm] = log_cholesky_estimation(acc,twist_avg,w,t)
    %theta_init = [-log(2)/2 0 0 0 0 0 0 0 0 0]';
    theta_init = 0.1 * randn(10,1); % small random around zero
    theta_init(1) = -log(2)/2 + 0.1 * randn(1); % keep alpha around -log(2)/2
    
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxIterations',50);
    %options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIterations',50);
    [theta, res_norm, b, c, output]  = lsqnonlin(@(theta) getResidual(acc,twist_avg,w,t,func_pi(theta)),theta_init,[],[],options);
    iterations = output.iterations;
    
    pi_ = func_pi(theta);

% [t2, twist_pred] = ode23t(@(t, twist) forwardDynamics(twist,t,w,piToInertiaMatrix(pi_)),timeSpan,initialTwist);
% plot(t,twist_avg, Marker=".")
% input("show pred")
% hold on
% plot(t2,twist_pred)
% legend('twist actual 1','twist actual 2','twist actual 3','twist actual 4','twist actual 5','twist actual 6',...
%        'twist pred 1','twist pred 2','twist pred 3','twist pred 4','twist pred 5','twist pred 6')
% xlabel('t');
% ylabel('twist');
% 
% input("continue")
% 
% clf;
end


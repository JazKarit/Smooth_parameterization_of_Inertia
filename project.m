clear;clc;clear all;


residuals = [];
iterations = [];

for i=1:20
    pi_actual = zeros(10,1);
    pi_actual(1) = sample_rand(0.01,1);
    pi_actual(2) = sample_rand(-.5,.5);
    pi_actual(3) = sample_rand(-.5,.5);
    pi_actual(4) = sample_rand(-.5,.5);
    
    L = [sample_rand(0.01,1);
         sample_rand(0.01,1);
         sample_rand(0.01,1)];
    
    d = [L(2)+L(3), L(1)+L(3), L(1)+L(2)]';
    
    D = diag(d);
    
    
    
    
    beta = [rand(1),rand(1),rand(1)];
    beta = beta / norm(beta);
    theta_r = rand(1)^(1/3) * 2 * pi;
    R = axisangle2rot(beta,theta_r);
    
    I_rot = R*D*R';
    
    pi_actual(5) = I_rot(1,1);
    pi_actual(6) = I_rot(2,2);
    pi_actual(7) = I_rot(3,3);
    pi_actual(8) = I_rot(1,2);
    pi_actual(9) = I_rot(2,3);
    pi_actual(10) = I_rot(1,3);
    
    h = pi_actual(2:4)';
    I = [I_rot, skewSymmetric(h);
         skewSymmetric(h)', pi_actual(1)*eye(3)];
    
    w = @(t) [sin(t); sin(t + pi/3); sin(t + 2*pi/3); sin(t + 3*pi/3); sin(t + 4*pi/3); sin(t + 5*pi/3)];
    
    timeSpan = [0 2*pi];
    initialTwist = [0;0;0;0;0;0];
    [t, twist] = ode23t(@(t, twist) forwardDynamics(twist,t,w,I),timeSpan,initialTwist);
    
    
    %noise = sample_rand(.075,.125,size(twist));
    
    noise = sample_rand(-0.001,0.001,size(twist));
    
    twist_noisy = twist;% + noise;
    twist_noisy = twist_noisy';
    twist = twist';
    
    
    
    acc = zeros(6,size(twist_noisy,2)-1);
    twist_avg = zeros(6,size(twist_noisy,2)-1);
    
    for i = 1:size(twist_noisy,2)-1
        acc(:,i) = (twist_noisy(:,i+1) - twist_noisy(:,i)) / (t(i+1)-t(i));
        twist_avg(:,i) = (twist_noisy(:,i+1) + twist_noisy(:,i)) / 2;
    end
    
    pi_init = [1 0 0 0 1 1 1 0 0 0]';
    
    theta_init = [-log(2)/2 0 0 0 0 0 0 0 0 0]';
    %func_pi(theta_init)
    
    
    
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIterations',50);
    %options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxIterations',50);
    %[pi_, a, b, c, output]  = lsqnonlin(@(pi) getResidual(acc,twist_avg,w,t,pi_),pi_init,[],[],options);
    [theta, a, b, c, output]  = lsqnonlin(@(theta) getResidual(acc,twist_avg,w,t,func_pi(theta)),theta_init,[],[],options);
    output.iterations;
    
    pi_ = func_pi(theta);
    
    residuals = [residuals, pi_actual-pi_];
    iterations = [iterations, output.iterations];
    

% [t2, twist_pred] = ode23t(@(t, twist) forwardDynamics(twist,t,w,piToInertiaMatrix(pi_)),timeSpan,initialTwist);
% plot(t,twist, Marker=".")
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

histogram(iterations)


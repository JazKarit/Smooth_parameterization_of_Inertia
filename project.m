%% 

method = 'log-cholesky';
method = 'normal';

%%
clear;clc;clear all;


pi_err = [];
iterations = [];
res_norms = [];

for i=1:1000
    pi_actual = generateRandomParameters();
    
    h = pi_actual(2:4)';

    I_rot = [pi_actual(5) pi_actual(8) pi_actual(10);
             pi_actual(8) pi_actual(6) pi_actual(9);
             pi_actual(10) pi_actual(9) pi_actual(7);];

    I = [I_rot, skewSymmetric(h);
         skewSymmetric(h)', pi_actual(1)*eye(3)];
    
    w = @(t) [sin(t); sin(t + pi/3); sin(t + 2*pi/3); sin(t + 3*pi/3); sin(t + 4*pi/3); sin(t + 5*pi/3)];
    
    timeSpan = [0 2*pi];
    initialTwist = [0;0;0;0;0;0];

    [t, twist] = ode23t(@(t, twist) forwardDynamics(twist,t,w,I),timeSpan,initialTwist);
    
    
    noise = sample_rand(.075,.125,size(twist));
    
    %noise = sample_rand(-0.001,0.001,size(twist));
    
    twist_noisy = twist + noise;
    twist_noisy = twist_noisy';
    twist = twist';
    
    % Define a uniform time grid
    dt = 0.1; % or any fixed step you want (matching the paper)
    t_uniform = (0:dt:2*pi)'; % column vector
    
    % Interpolate twist onto the uniform grid
    twist_noisy_interp = interp1(t, twist_noisy', t_uniform, 'linear')'; % interpolate each row
    twist_interp = interp1(t, twist', t_uniform, 'linear')'; % clean (optional)
    
    acc = diff(twist_noisy,1,2) / dt;
    twist_avg = (twist_noisy(:,1:end-1) + twist_noisy(:,2:end))/2;

    
    [pi_, iterations_,res_norm] = log_cholesky_estimation(acc,twist_avg,w,t);

    %[pi_, iterations_,res_norm] = regular_estimation(acc,twist_avg,w,t);

    %[pi_, iterations_, res_norm] = eigenvalue_estimation(acc,twist_avg,w,t);
    
    %pi_actual-pi_
    pi_err = [pi_err, norm(pi_actual-pi_)];
    res_norms = [res_norms, res_norm];

    % if norm(pi_actual-pi_) > 0.1
    %     iterations_ = 52;
    % end

    iterations = [iterations, iterations_];
    
   

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
% histogram(pi_err)



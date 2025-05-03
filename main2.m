function RunInertiaEstimation()
    % Compare inertial parameter estimation methods on simulated data
    % Implementation of "Smooth Parameterization of Rigid-Body Inertia"
    
    clear; clc;
    
    % Storage for results
    all_iterations = cell(3,1);
    for i = 1:3
        all_iterations{i} = [];
    end
    
    % Run 1000 simulations as in the paper
    num_sims = 100; % Reduced for testing - use 1000 for full replication
    for sim = 1:num_sims
        fprintf('Simulation %d of %d\n', sim, num_sims);
        
        % Generate random physically consistent inertia parameters
        pi_true = generateRandomParameters();
        
        % Get the corresponding spatial inertia tensor
        I_true = piToInertiaMatrix(pi_true);
        
        % Generate simulation data
        [t, twist, twist_noisy, acc, twist_midpoint, wrench] = simulateRigidBody(I_true);
        
        % Method 1: Log-Cholesky Parameterization
        theta_init = zeros(10,1);
        options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', ...
                               'MaxIterations', 50, 'Display', 'off');
        [~, ~, ~, ~, output1] = lsqnonlin(@(theta) residualLogCholesky(theta, acc, twist_midpoint, wrench), ...
                                         theta_init, [], [], options);
        all_iterations{1} = [all_iterations{1}, output1.iterations];
        
        % Method 2: Eigenvalue Decomposition Parameterization
        phi_init = zeros(10,1);
        phi_init(1) = 1; % mass = 1
        phi_init(5:7) = 0.1; % L = [0.1, 0.1, 0.1]
        [~, ~, ~, ~, output2] = lsqnonlin(@(phi) residualEigenvalue(phi, acc, twist_midpoint, wrench), ...
                                         phi_init, [], [], options);
        all_iterations{2} = [all_iterations{2}, output2.iterations];
        
        % Method 3: Exponential Map Parameterization
        Q_init = zeros(10,1);
        [~, ~, ~, ~, output3] = lsqnonlin(@(Q) residualExponential(Q, acc, twist_midpoint, wrench), ...
                                         Q_init, [], [], options);
        all_iterations{3} = [all_iterations{3}, output3.iterations];
    end
    
    % Plot histograms similar to Fig. 1 in the paper
    methods = {'Log-Cholesky', 'Eigenvalue', 'Exponential'};
    figure('Position', [100, 100, 800, 600]);
    
    for i = 1:3
        subplot(3,1,i);
        histogram(all_iterations{i}, 'BinWidth', 1);
        title(methods{i});
        xlabel('Number of Iterations');
        ylabel('Count');
        xlim([0, 51]);
    end
    
    sgtitle('Convergence Performance of Inertia Parameterizations');
end

function pi_true = generateRandomParameters()
    % Generate random physically consistent inertial parameters
    % as described in Table I of the paper
    
    % Mass: uniform in [0.01, 1]
    m = rand() * 0.99 + 0.01;
    
    % Center of mass: uniform in [-0.5, 0.5] for each component
    c = (rand(3,1) - 0.5);
    
    % Second moments: uniform in [0.01, 1] for each
    L = rand(3,1) * 0.99 + 0.01;
    
    % Convert second moments to principal moments using equation (5)
    D = [L(2)+L(3); L(1)+L(3); L(1)+L(2)];
    
    % Random rotation
    % Generate random axis and angle
    axis = randn(3,1);
    axis = axis / norm(axis);
    angle = rand() * pi; % Random angle between 0 and π
    
    % Convert to rotation matrix using Rodrigues' formula
    R = axisangle2rot(axis, angle);
    
    % Compute rotational inertia tensor
    I_rot = R * diag(D) * R';
    
    % Construct true parameter vector
    h = m * c;
    pi_true = [m; h; I_rot(1,1); I_rot(2,2); I_rot(3,3); I_rot(1,2); I_rot(2,3); I_rot(1,3)];
end

function [t, twist, twist_noisy, acc, twist_midpoint, wrench] = simulateRigidBody(I)
    % Simulate rigid body dynamics with the given inertia matrix
    
    % Define wrench function as in the paper
    w = @(t) [sin(t); sin(t + pi/3); sin(t + 2*pi/3); 
              sin(t + pi); sin(t + 4*pi/3); sin(t + 5*pi/3)];
    
    % Simulation time span
    timeSpan = [0, 2*pi];
    
    % Initial state (zero twist)
    initialTwist = zeros(6,1);
    
    % Simulate rigid body dynamics
    [t, twist_raw] = ode45(@(t, twist) forwardDynamics(twist, t, w, I), timeSpan, initialTwist);
    
    % Convert to matrix form
    twist = twist_raw';
    
    % Resample to uniform time step of 0.1s
    dt = 0.1;
    t_uniform = (0:dt:2*pi)';
    twist_uniform = interp1(t, twist', t_uniform, 'linear')';
    
    % Generate noisy measurements with bias as described in the paper
    % Uniform distribution with width 0.05 and mean 0.1
    noise = (rand(size(twist_uniform)) - 0.5) * 0.05 + 0.1;
    twist_noisy = twist_uniform + noise;
    
    % Compute accelerations using midpoint method
    acc = diff(twist_noisy, 1, 2) / dt;
    
    % Compute midpoint twist
    twist_midpoint = (twist_noisy(:,1:end-1) + twist_noisy(:,2:end)) / 2;
    
    % Corresponding time points
    t = t_uniform(1:end-1) + dt/2;
    
    % Compute wrench at midpoint times
    wrench = zeros(6, length(t));
    for i = 1:length(t)
        wrench(:,i) = w(t(i));
    end
    
    % Return only the uniform time samples
    t = t_uniform;
end

function dtwist = forwardDynamics(twist, t, wrenchFcn, I)
    % Forward dynamics for a rigid body
    % twist: [ω; v] 6x1 spatial velocity
    % I: 6x6 spatial inertia matrix
    % wrenchFcn: function handle that returns wrench at time t
    
    % Extract angular and linear velocity
    omega = twist(1:3);
    
    % Compute cross product operator
    crossOp = [skewSymmetric(omega), skewSymmetric(twist(4:6));
               zeros(3), skewSymmetric(omega)];
    
    % Compute current wrench
    wrench = wrenchFcn(t);
    
    % Compute acceleration
    dtwist = I \ (wrench - crossOp * I * twist);
end

function S = skewSymmetric(v)
    % Convert vector to skew-symmetric matrix
    S = [0, -v(3), v(2); 
         v(3), 0, -v(1); 
         -v(2), v(1), 0];
end

function R = axisangle2rot(axis, angle)
    % Convert axis-angle representation to rotation matrix
    axis = axis / norm(axis); % Ensure unit vector
    
    % Rodrigues' rotation formula
    K = skewSymmetric(axis);
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;
end

function I = piToInertiaMatrix(pi)
    % Convert inertial parameter vector to spatial inertia matrix
    
    m = pi(1);
    h = pi(2:4);
    Ixx = pi(5); Iyy = pi(6); Izz = pi(7);
    Ixy = pi(8); Iyz = pi(9); Ixz = pi(10);
    
    I_rot = [Ixx, Ixy, Ixz;
             Ixy, Iyy, Iyz;
             Ixz, Iyz, Izz];
         
    % Construct spatial inertia matrix
    I = [I_rot, skewSymmetric(h);
         skewSymmetric(h)', m*eye(3)];
end

function res = residualLogCholesky(theta, acc, twist, wrench)
    % Residual function for log-Cholesky parameterization
    
    % Convert parameters to inertial parameters using equation (9)
    pi = logCholeskyToPi(theta);
    
    % Compute spatial inertia
    I = piToInertiaMatrix(pi);
    
    % Compute residuals for each time point
    res = zeros(size(wrench));
    for i = 1:size(twist, 2)
        xi = twist(:,i);
        xi_dot = acc(:,i);
        
        % Compute cross operator
        omega = xi(1:3);
        v = xi(4:6);
        xi_cross = [skewSymmetric(omega), skewSymmetric(v);
                    zeros(3), skewSymmetric(omega)];
        
        % Equation (13) from the paper
        res(:,i) = I * xi_dot + xi_cross * I * xi - wrench(:,i);
    end
    
    % Flatten residual vector
    res = res(:);
end

function pi = logCholeskyToPi(theta)
    % Convert log-Cholesky parameters to inertial parameters
    % Implementation of equation (9) from the paper
    
    alpha = theta(1);
    d = theta(2:4);
    s = [theta(5), theta(7), theta(6)]; % [s12, s13, s23]
    t = theta(8:10);
    
    % Compute exponentials
    e_alpha = exp(alpha);
    e_d = exp(d);
    
    % Construct parameters (equation 9)
    pi = zeros(10,1);
    
    % Mass
    pi(1) = e_2alpha * (t(1)^2 + t(2)^2 + t(3)^2 + 1);
    
    % First mass moments
    pi(2) = e_2alpha * t(1) * e_d(1);
    pi(3) = e_2alpha * (t(1)*s(1) + t(2)*e_d(2));
    pi(4) = e_2alpha * (t(1)*s(2) + t(2)*s(3) + t(3)*e_d(3));
    
    % Diagonal elements of inertia tensor
    pi(5) = e_2alpha * (s(1)^2 + s(2)^2 + s(3)^2 + e_d(2)^2 + e_d(3)^2);
    pi(6) = e_2alpha * (s(2)^2 + s(3)^2 + e_d(1)^2 + e_d(3)^2);
    pi(7) = e_2alpha * (s(1)^2 + e_d(1)^2 + e_d(2)^2);
    
    % Off-diagonal elements of inertia tensor
    pi(8) = -e_2alpha * s(1) * e_d(1);
    pi(9) = -e_2alpha * (s(1)*s(2) + s(3)*e_d(2));
    pi(10) = -e_2alpha * s(2) * e_d(1);
    
    function val = e_2alpha
        val = exp(2*alpha);
    end
end

function res = residualEigenvalue(phi, acc, twist, wrench)
    % Residual function for eigenvalue parameterization
    
    % Extract parameters
    m = phi(1);
    c = phi(2:4);
    beta = phi(5:7); % Rotation parameters
    L = phi(8:10);   % Second moments
    
    % Make sure second moments are positive
    L = abs(L);
    
    % Construct principal moments
    D = [L(2)+L(3); L(1)+L(3); L(1)+L(2)];
    
    % Calculate rotation matrix from exponential coordinates
    beta_norm = norm(beta);
    if beta_norm < 1e-10
        R = eye(3);
    else
        axis = beta / beta_norm;
        angle = beta_norm;
        R = axisangle2rot(axis, angle);
    end
    
    % Calculate rotational inertia
    I_rot = R * diag(D) * R';
    
    % First mass moment
    h = m * c;
    
    % Construct inertial parameter vector
    pi = [m; h; I_rot(1,1); I_rot(2,2); I_rot(3,3); I_rot(1,2); I_rot(2,3); I_rot(1,3)];
    
    % Compute spatial inertia
    I = piToInertiaMatrix(pi);
    
    % Compute residuals
    res = zeros(size(wrench));
    for i = 1:size(twist, 2)
        xi = twist(:,i);
        xi_dot = acc(:,i);
        
        % Compute cross operator
        omega = xi(1:3);
        v = xi(4:6);
        xi_cross = [skewSymmetric(omega), skewSymmetric(v);
                    zeros(3), skewSymmetric(omega)];
        
        % Equation (13) from the paper
        res(:,i) = I * xi_dot + xi_cross * I * xi - wrench(:,i);
    end
    
    % Flatten residual vector
    res = res(:);
end

function res = residualExponential(Q, acc, twist, wrench)
    % Residual function for exponential map parameterization
    
    % Convert parameters to a symmetric 4x4 matrix
    S4 = zeros(4);
    idx = 1;
    for i = 1:4
        for j = i:4
            if i == j
                S4(i,j) = Q(idx);
            else
                S4(i,j) = Q(idx) / sqrt(2);
                S4(j,i) = Q(idx) / sqrt(2);
            end
            idx = idx + 1;
            if idx > 10
                break;
            end
        end
    end
    
    % Calculate pseudo-inertia using matrix exponential
    J = expm(S4);
    
    % Extract inertial parameters from J
    m = J(4,4);
    h = J(1:3,4);
    Sigma = J(1:3,1:3);
    I_rot = trace(Sigma)*eye(3) - Sigma;
    
    % Construct parameter vector
    pi = [m; h; I_rot(1,1); I_rot(2,2); I_rot(3,3); I_rot(1,2); I_rot(2,3); I_rot(1,3)];
    
    % Compute spatial inertia
    I = piToInertiaMatrix(pi);
    
    % Compute residuals
    res = zeros(size(wrench));
    for i = 1:size(twist, 2)
        xi = twist(:,i);
        xi_dot = acc(:,i);
        
        % Compute cross operator
        omega = xi(1:3);
        v = xi(4:6);
        xi_cross = [skewSymmetric(omega), skewSymmetric(v);
                    zeros(3), skewSymmetric(omega)];
        
        % Equation (13) from the paper
        res(:,i) = I * xi_dot + xi_cross * I * xi - wrench(:,i);
    end
    
    % Flatten residual vector
    res = res(:);
end
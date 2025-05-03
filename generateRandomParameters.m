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
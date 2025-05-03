function pi = phi_to_pi(phi)
%FUNC_PHI Summary of this function goes here
%   Detailed explanation goes here
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
end


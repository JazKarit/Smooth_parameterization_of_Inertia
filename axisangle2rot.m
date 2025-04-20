function R = axisangle2rot(omega,theta)
 % axisangle2rot: Given the axis (omega) and angle (theta) 
 %                of rotation, returns the corresponding rotation matrix
 %                theta is a scalar 
 %                omega is a unit vector
     W = skewSymmetric(omega);
     R = eye(3) + sin(theta)*W + (1-cos(theta))*W^2 ;
end
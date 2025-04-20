function I = piToInertiaMatrix(pi)
%piITOINERTIAMATRIX Summary of this function goes here
%   Detailed expilanation goes here
I_rot = [pi(5) pi(8) pi(10);
         pi(8) pi(6) pi(9);
         pi(10) pi(9) pi(7)];

h = pi(2:4)';

I = [I_rot, skewSymmetric(h);
     skewSymmetric(h)', pi(1)*eye(3)];

end


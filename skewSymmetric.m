function W = skewSymmetric(w)
% skewSymmetric: Given a vector in R^3 returns the corresponding
%                skew symmetric matrix (3x3)
        W = [0, -w(3), w(2);
             w(3), 0, -w(1);
             -w(2), w(1), 0];
 end
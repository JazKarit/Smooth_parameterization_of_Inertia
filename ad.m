function adV = ad(V)
% Given a twist V1, calculate the corresponding 6x6 [adV1] matrix.
% This matrix can be used to calculate the Lie Bracket of V1 and another twist V2
% by ad(V1)*V2
    w = V(1:3);
    v = V(4:6);
    adV = [skewSymmetric(w) zeros(3,3);
           skewSymmetric(v) skewSymmetric(w)];
end

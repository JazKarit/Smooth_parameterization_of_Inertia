function residual = getResidual(acc, twist, w, t, pi)
%GETRESIDUAL Summary of this function goes here
%   Detailed explanation goes here

I = piToInertiaMatrix(pi);

w_pred = zeros(6,size(twist,2));
residual = zeros(6,size(twist,2));

for i = 1:size(twist,2)
    w_pred(:,i) = I*acc(:,i)+ad(twist(:,i))'*I*twist(:,i);
    residual(:,i) = w_pred(:,i)-w((t(i+1)+t(i))/2); % t(i) or t(i+0.5)
end
reshape(residual.',1,[]);

end


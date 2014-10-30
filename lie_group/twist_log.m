function [W] = twist_log(M)

R = M(1:3, 1:3);
val = 0.5*(trace(R)-1);

% This is 
if val > 1.0
    val = 1.0;
elseif val < -1.0
    val =-1.0;
end

theta = acos(val);

if 0==theta
    w_hat = zeros(3,3);
else
    w_hat = (R-R')/(2*sin(theta))*theta;
end

if theta < 1e-1
    V_inv = eye(3) - 0.5*w_hat + 1/12*w_hat*w_hat;
else
    V_inv = eye(3) - 0.5*w_hat + (1-0.5*theta*sin(theta)/(1-cos(theta)))/theta^2*w_hat*w_hat;
end

u = V_inv*M(1:3, 4);

W = [so3_gla(w_hat); u];
end

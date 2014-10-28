function [w] = screw_log(R)

val = 0.5*(trace(R)-1);

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

w = so3_gla(w_hat);
end
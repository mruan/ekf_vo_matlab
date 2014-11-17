function [w] = screw_log(R)

val = 0.5*(trace(R)-1);

if val > 1.0
    val = 1.0;
elseif val < -1.0
    val =-1.0;
end

theta = acos(val);

if 0==theta
    w = zeros(3,1);
elseif pi==theta
    K2 = 0.5*(R+eye(3));
    w = sqrt(diag(K2));
    if     K2(1,1) ~= 0 % if w1 ~=0
        w = theta*w.*[      1      ; sign(K2(1,2)); sign(K2(1,3))];
    elseif K2(2,2) ~= 0 % if w2 ~=0
        w = theta*w.*[sign(K2(2,1));      1       ; sign(K2(2,3))];
    else                % if w3 ~=0
        w = theta*w.*[sign(K2(3,1)); sign(K2(3,2));      1       ];
    end
else
    w_hat = (R-R')/(2*sin(theta))*theta;
    w = so3_gla(w_hat);
end

end
function [z_hn] = proj_P3(x, K, R, t)
%PROJ_P3 Summary of this function goes here
%   Detailed explanation goes here

z_h = K*(R*x(1:3)+t*x(4));
            
% z_hn i.e. \hat(z) (with normalization in image coord. now)
z_hn = z_h/z_h(3);
z_hn = z_hn(1:2);
end


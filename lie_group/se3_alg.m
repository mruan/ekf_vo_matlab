function [ alg ] = se3_alg( w )
%ALG_SE3 Summary of this function goes here
%   Detailed explanation goes here

assert(numel(w) == 6, 'w must have exactly 6 elements');
% Assume the first 3 elements are linear speed, then come the angular vel
wx = so3_alg(w(4:6));

alg = [wx w(1:3); 0 0 0 0];
end


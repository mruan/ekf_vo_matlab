function [W] = se3_gla(M)
% inverse of alg operator for se3 (alg -> gla)

% Assert the last row of M is [0 0 0 0]
assert(norm(M(4,:)) == 0, 'Last row of M is not [0 0 0 0]');

% Assume the first 3 elements are linear speed, then come the angular vel
W = [M(1:3, 4); so3_gla(M(1:3, 1:3))];

end
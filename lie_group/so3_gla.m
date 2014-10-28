function [vec] = so3_gla(m)
% inverse of alg operator for so3 (alg -> gla)

assert(norm(m+transpose(m)) == 0, 'm is not a symmetric-skew matrix');

vec = [m(3,2); m(1,3); m(2,1)];
end
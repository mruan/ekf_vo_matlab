function [vec] = so3_gla(m)
% inverse of alg operator for so3 (alg -> gla)

assert(norm(m+transpose(m)) < 1e-10, 'm is not a symmetric-skew matrix');

vec = 0.5*[m(3,2)-m(2,3); m(1,3)-m(3,1); m(2,1)-m(1,2)];
end
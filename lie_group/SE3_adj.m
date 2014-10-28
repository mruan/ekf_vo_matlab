function [Adj] = SE3_adj(M)

R = M(1:3, 1:3);
t = M(1:3, 4);

% Adj = [R  t_x R; 0 R]
Adj = [R so3_alg(t)*R; zeros(3) R];

end
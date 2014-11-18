function [Mx] = so3_alg(v)

assert(numel(v) == 3, 'Input vector must contain exactly 3 elements');

Mx = [   0 -v(3)  v(2);...
      v(3)     0 -v(1);...
     -v(2)  v(1)    0];
end
% Test the following derivative:
%   d/dx exp(x)'*a
%  =d/dx exp(-x)*a
%  =d/dx exp(-(z<+>x))
%  =d/dz exp(-x)*exp(-z)*a
% ...
%  = [exp(x)'*a]_x*exp(x)'

rng(1);
%
rx = rand(3,1);
a = rand(3,1);
Rx = screw_exp(rx);

delta = 1e-6;
Hn = zeros(3,3);

for i=1:3
   d_rx = [0 0 0]';
   d_rx(i) = delta;
   
   Rx_p = screw_exp(d_rx)*Rx;
   
   ai = Rx'  *a;
   ao = Rx_p'*a;
   
   Hn(:,i) = (ao-ai)/delta;
end

Hn

so3_alg(Rx'*a)*Rx'

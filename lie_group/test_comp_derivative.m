

rng(0);

K = [500 0 320; 0 500 240; 0 0 1];

x = 10*rand(3,1);

theta_1 = rand(3,1);

r1 = eye(3);%screw_exp(theta_1);

T_x = K*r1*x;

delta = 1e-6;
for i=1:3
theta = zeros(3,1);
theta(i) = delta;
T_dx = K*screw_exp(theta)*r1*x;

dx = T_dx - T_x;

H(:,i) = dx /delta;
end

G1 = [0 0 0; 0 0 -1; 0 1 0];
G2 = [0 0 1; 0 0 0; -1 0 0];
G3 = [0 -1 0; 1 0 0; 0 0 0];

H
Hp = [K*G1*r1*x K*G2*r1*x K*G3*r1*x]
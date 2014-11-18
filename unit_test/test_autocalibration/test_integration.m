% Test Numerical Integration of a trajectory path
% serves as a unit test to test_autocalib.m
% test_sum is a unit test to this program

%{
a = sym('a', [1 1]);
a_dot = sym('a_dot', [1 1]);
a_ddot = sym('a_ddot', [1 1]);

R=[
[1     0       0 ]
[0 cos(a) -sin(a)]
[0 sin(a)  cos(a)]];

T=[0; cos(a); sin(a)]*r
T_dot = [0; -sin(a); cos(a)]*r*a_dot
(Just for a comparison: 
T_ddot = [0; -cos(a); -sin(a)]*r*a_dot^2 + [0; -sin(a); cos(a)]*r*a_ddot
}
R_dot=[
[0      0       0 ]
[0 -sin(a) -cos(a)]
[0  cos(a) -sin(a)]]*a_dot;

w_sb_b_x = R'*R_dot
w_sb_b_x = [
[ 0  0  0]
[ 0  0 -1]
[ 0  1  0]]*a_dot

v_sb_b = R'*T_dot
v_sb_b = [0 0 1]'*r*a_dot

v_sb_b_dot = [0 0 1]'*r*a_ddot

a_sb_b = v_sb_b_dot + w_sb_b_x * v_sb_b;
a_sb_b = [0; -r*a_dot^2; r*a_ddot];
%}
function test_integration
dt = 0.01;
t  = 0:dt:10;
r  = 1.0;

% a      = t;
% a_dot  =  ones(size(t));
% a_ddot = zeros(size(t));
% 
% a      = 0.5*t.^2;
% a_dot  =     t;
% a_ddot = ones(size(t));

a      = 0.01*t.^3 + 0.05*t.^2 + t;
a_dot  = 0.03*t.^2 +  0.1*t    + 1;
a_ddot = 0.06*t    +  0.1;

N = numel(t);

rx_gt = zeros(3, N);
tx_gt = zeros(3, N);
rx_it = zeros(3, N);
tx_it = zeros(3, N);
for i=1:N
   R = [1 0 0; 0 cos(a(i)) -sin(a(i)); 0 sin(a(i)) cos(a(i))];
   
   rx_gt(:,i) = screw_log(R);
   tx_gt(:,i) = [0; cos(a(i)); sin(a(i))]*r;
end

rx_it(:,1) = rx_gt(:,1);
tx_it(:,1) = tx_gt(:,1);
R_im = screw_exp(rx_it(:, 1));
vt = [0; -sin(a(1)); cos(a(1))]*r*a_dot(1);
at = R_im*[0; -r*a_dot(1)^2; r*a_ddot(1)^2];
for i=2:N
    wx = [1 0 0]'*a_dot(i-1);
%     ux = [1 0 0]'*a_ddot(i-1);
    
    rx_it(:,i) = screw_log(R_im*screw_exp(wx*dt));
    tx_it(:,i) = tx_it(:,i-1) + vt*dt;
%     vt = [0; -sin(a(i)); cos(a(i))]*r*a_dot(i);
    vt = vt + at*dt;
    
    R_im = screw_exp(rx_it(:,i));
    at = R_im*[0; -r*a_dot(i)^2; r*a_ddot(i)];
end

figure(1);
subplot(2,1,1);
plot(t, rx_gt(1,:), 'r',   t, rx_gt(2,:), 'g',   t, rx_gt(3,:), 'b'); hold on;
plot(t, rx_it(1,:), 'r-.', t, rx_it(2,:), 'g-.', t, rx_it(3,:), 'b-.'); hold off;
subplot(2,1,2);
plot(t, tx_gt(1,:), 'r',   t, tx_gt(2,:), 'g',   t, tx_gt(3,:), 'b'); hold on;
plot(t, tx_it(1,:), 'r-.', t, tx_it(2,:), 'g-.', t, tx_it(3,:), 'b-.'); hold off;

[max_rx, id]  = max(abs(rx_gt(1,:) - rx_it(1,:)));
end
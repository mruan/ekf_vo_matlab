% Test Numerical Integration of a scalar
% As a unit test to test_integration

dt = 0.01;
t  = 0:dt:4*pi;
r  = 1.0;

a      = 0.1*t.^3 + 0.5*t.^2 + 1;
a_dot  = 0.3*t.^2 +     t;
a_ddot = 0.6*t + 1;

N = numel(t);

A = zeros(size(a));
A_dot =zeros(size(A));
% A_ddot= zeros(size(A));

A(1) = a(1);
A_dot(1) = a_dot(1);

for i=2:N
    A(i) = A(i-1) + A_dot(i-1)*dt;
    A_dot(i) = A_dot(i-1) + a_ddot(i-1)*dt;
end

subplot(2,1,1)
plot(t, a, 'r', t, A, 'b');
subplot(2,1,2)
plot(t, a_dot, 'r', t, A_dot, 'b');
% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

clear;

load('sim_env.mat'); % measurement

t = measurement.Time;
N = numel(t);

rx = zeros(3, N); rx(:,1) = measurement.True(1, 1:3)';
tx = zeros(3, N); tx(:,1) = measurement.True(1, 4:6)';
vt = zeros(3, N); vt(:,1) = measurement.True(1, 10:12)';

v = vt(:,1);
w = measurement.True(1, 7:9)';
a = measurement.True(1, 13:15)';

R = screw_exp(rx(:,1));
for i=2:N
    
    dt = t(i) - t(i-1);
    R = R*screw_exp(w*dt);
    
    rx(:,i) = screw_log(R);
    
    tx(:,i) = tx(:,i-1) + v*dt;
    
    v = v + a*dt;
    
    vt(:,i) = v;
    
    a = measurement.True(i, 13:15)';
    w = measurement.True(i, 7:9)';
end

subplot(3,1,1);
plot(t, rx(1,:)'-measurement.True(:, 1), 'r', ...
     t, rx(2,:)'-measurement.True(:, 2), 'g', ...
     t, rx(3,:)'-measurement.True(:, 3), 'b');
 
subplot(3,1,2);
plot(t, tx(1,:)'-measurement.True(:, 4), 'r', ...
     t, tx(2,:)'-measurement.True(:, 5), 'g', ...
     t, tx(3,:)'-measurement.True(:, 6), 'b');
 
subplot(3,1,3);
plot(t, vt(1,:)'-measurement.True(:, 10), 'r', ...
     t, vt(2,:)'-measurement.True(:, 11), 'g', ...
     t, vt(3,:)'-measurement.True(:, 12), 'b');
 
% subplot(3,1,1);
% plot(t, measurement.True(:, 1), 'r', ...
%      t, measurement.True(:, 2), 'g', ...
%      t, measurement.True(:, 3), 'b'); hold on;
% plot(t, rx(1,:), 'r-.', t, rx(2,:), 'g-.', t, rx(3,:), 'b-.'); hold off;
% 
% subplot(3,1,2);
% plot(t, measurement.True(:, 4), 'r', ...
%      t, measurement.True(:, 5), 'g', ...
%      t, measurement.True(:, 6), 'b'); hold on;
% plot(t, tx(1,:), 'r-.', t, tx(2,:), 'g-.', t, tx(3,:), 'b-.'); hold off;
% 
% subplot(3,1,3);
% plot(t, measurement.True(:, 10), 'r', ...
%      t, measurement.True(:, 11), 'g', ...
%      t, measurement.True(:, 12), 'b'); hold on;
% plot(t, vt(1,:), 'r-.', t, vt(2,:), 'g-.', t, vt(3,:), 'b-.'); hold off;

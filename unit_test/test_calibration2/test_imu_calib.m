% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

%clear;
function [t, f, rx, tx, vt, at] = test_imu_calib()
load('sim_env.mat'); % measurement

t = measurement.Time;
N = numel(t);

rx = zeros(3, N);
tx = zeros(3, N);
vt = zeros(3, N);
at = zeros(3, N);

f = ekf_autocalib();
f.rsb = measurement.True(1, 1:3)';
f.tsb = measurement.True(1, 4:6)';
%f.wt  = measurement.True(1, 7:9)';
%f.vt  = measurement.True(1, 10:12)';

R = screw_exp(rx(:,1));
for i=1:N
    w = measurement.Data(i, 1:3)';
    a = measurement.Data(i, 4:6)';
    
    imu_noise = 1e-6*eye(6);
    
    f.onImuUpdate(w, a, imu_noise, t(i));
    
    rx(:,1) = f.rsb;
    tx(:,1) = f.tsb;
    vt(:,1) = f.vt;
    wt(:,1) = f.wt;
    at(:,i) = f.at;
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
end

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

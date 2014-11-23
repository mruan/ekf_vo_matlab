% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

%clear;
function [f] = test_imu_calib()
dbstop if error;

load('sim_env3.mat'); % measurement
addpath('../../lie_group');

t = measurement.Time;
N = numel(t);

f = ekf_imucalib();
X = zeros(f.N_states,1);
P = 100*eye(f.N_states);
log_X = zeros(f.N_states, N);
log_P = zeros(f.N_states, N);

% f.X(f.rsb) = measurement.True(1, 1:3)'; f.P(f.rsb, f.rsb) = 1e-0*eye(3);
% f.X(f.tsb) = measurement.True(1, 4:6)'; f.P(f.tsb, f.tsb) = 1e-0*eye(3);
% f.X(f.wt)  = measurement.True(1, 7:9)'; f.P(f.wt, f.wt) = 1e-0*eye(3);
% % f.X(f.vt)  = measurement.True(1, 10:12)'; f.P(f.vt, f.vt) = 1e-8*eye(3);
% % f.X(f.at)  = measurement.True(1, 13:15)'; f.P(f.at, f.at) = 1e-8*eye(3);
% % f.X(f.g0)  = [0 0 9.8]'; f.P(f.g0, f.g0) = eye(3);
% f.X(f.wb) = [0.0 0.0 0.0]'; f.P(f.wb, f.wb) = 9e-2*eye(3);
% f.X(f.ab) = [0.0 0.0 0.0]'; f.P(f.ab, f.ab) = 9e-2*eye(3);
t_prev= 0;
for i=1:N
    t_now = t(i);
    [X_h, P_h] = f.propagate(X, P, t_now-t_prev);
    x_ref = measurement.True(i,:)';
    
    if measurement.Type(i) == 0
        w = measurement.Data(i, 1:3)';
        a = measurement.Data(i, 4:6)';
        imu_noise = 1e-6*eye(6);

        [X, P] = f.onImuUpdate(w, a, imu_noise, X_h, P_h);
        
    else % Camera update       
        ri = measurement.Data(i, 1:3)';
        ti = measurement.Data(i, 4:6)';
        
        omn_noise = 1e-6*eye(6);
            
        [X, P] = f.onOmnUpdate(ri, ti, omn_noise, X_h, P_h);
    end
    t_prev = t_now;
    
    log_X(:,i) = f.states_add(X, -x_ref);
    log_P(:,i) = sqrt(diag(P));
end

idx = 100:N;
subplot(4,2,1); plot_one(t(idx), log_X(f.rsb,idx), log_P(f.rsb,idx), 'rx');
subplot(4,2,2); plot_one(t(idx), log_X(f.tsb,idx), log_P(f.tsb,idx), 'tx');
subplot(4,2,3); plot_one(t(idx), log_X(f.wt, idx), log_P(f.wt, idx), 'wt');
subplot(4,2,4); plot_one(t(idx), log_X(f.vt, idx), log_P(f.vt, idx), 'vt');
subplot(4,2,5); plot_one(t(idx), log_X(f.at, idx), log_P(f.at, idx), 'at');
subplot(4,2,6); plot_one(t(idx), log_X(f.wb, idx), log_P(f.wb, idx), 'wb');
subplot(4,2,7); plot_one(t(idx), log_X(f.ab, idx), log_P(f.ab, idx), 'ab');
end

function [] = plot_one(t, x, p, name)
plot(t, x(1,:), 'r',   t, x(2,:), 'g',   t, x(2,:), 'b',...
     t, p(1,:), 'r-.', t, p(2,:), 'g-.', t, p(3,:), 'b-.',...
     t,-p(1,:), 'r-.', t,-p(2,:), 'g-.', t,-p(3,:), 'b-.'); title(name);
 
end
function [x] = extract_sub(x)
   x = x([1:12 16:18 19:24]);
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

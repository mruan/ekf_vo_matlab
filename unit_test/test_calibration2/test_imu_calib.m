% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

%clear;
function [f] = test_imu_calib()
dbstop if error;

load('sim_env_imu.mat'); % measurement
addpath('../../lie_group');

t = measurement.Time;
N = numel(t);

f = ekf_imucalib();
X = zeros(f.N_states,1);
P = 100*eye(f.N_states);
log_X = zeros(f.N_states, N);
log_P = zeros(f.N_states, N);

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

idx = 1:N;
subplot(4,2,1); plot_one(t(idx), log_X(f.rsb,idx), log_P(f.rsb,idx), 'rx');
subplot(4,2,2); plot_one(t(idx), log_X(f.tsb,idx), log_P(f.tsb,idx), 'tx');
subplot(4,2,3); plot_one(t(idx), log_X(f.wt, idx), log_P(f.wt, idx), 'wt');
subplot(4,2,4); plot_one(t(idx), log_X(f.vt, idx), log_P(f.vt, idx), 'vt');
subplot(4,2,5); plot_one(t(idx), log_X(f.at, idx), log_P(f.at, idx), 'at');
subplot(4,2,6); plot_one(t(idx), log_X(f.wb, idx), log_P(f.wb, idx), 'wb');
subplot(4,2,7); plot_one(t(idx), log_X(f.ab, idx), log_P(f.ab, idx), 'ab');
subplot(4,2,8); plot_one(t(idx), log_X(f.g0, idx), log_P(f.g0, idx), 'g0');
end

function [] = plot_one(t, x, p, name)
plot(t, x(1,:), 'r',   t, x(2,:), 'g',   t, x(2,:), 'b',...
     t, p(1,:), 'r-.', t, p(2,:), 'g-.', t, p(3,:), 'b-.',...
     t,-p(1,:), 'r-.', t,-p(2,:), 'g-.', t,-p(3,:), 'b-.'); title(name);
 
end

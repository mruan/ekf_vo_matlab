% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

%clear;
function [f] = test_dual_calib()
dbstop if error;

load('sim_env_dual.mat'); % measurement
addpath('../../lie_group');

t = measurement.Time;
N = numel(t);

f = ekf_dualcalib();
X = zeros(f.N_states,1);
P = 100*eye(f.N_states);

%% Initialize the relative pose of the camera in body frame
X(f.rbc) = [pi; 0; 0]; P(f.rbc, f.rbc) = 0.0025*eye(3);
X(f.tbc) = [0; 0; 0];  P(f.tbc, f.tbc) = 0.0025*eye(3);

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
        rsc = measurement.Data(i, 1:3)';
        tsc = measurement.Data(i, 4:6)';
        
        cam_noise = 1e-6*eye(6);
            
        [X, P] = f.onCamUpdate(rsc, tsc, cam_noise, X_h, P_h);
    end
    t_prev = t_now;
    
    log_X(:,i) = f.states_add(X, -x_ref);
    log_P(:,i) = sqrt(diag(P));
end

idx = 1:N;
subplot(5,2,1); plot_one(t(idx), log_X(f.rsb,idx), log_P(f.rsb,idx), 'rsb');
subplot(5,2,2); plot_one(t(idx), log_X(f.tsb,idx), log_P(f.tsb,idx), 'tsb');
subplot(5,2,3); plot_one(t(idx), log_X(f.wt, idx), log_P(f.wt, idx), 'wt');
subplot(5,2,4); plot_one(t(idx), log_X(f.vt, idx), log_P(f.vt, idx), 'vt');
subplot(5,2,5); plot_one(t(idx), log_X(f.at, idx), log_P(f.at, idx), 'at');
subplot(5,2,6); plot_one(t(idx), log_X(f.wb, idx), log_P(f.wb, idx), 'wb');
subplot(5,2,7); plot_one(t(idx), log_X(f.ab, idx), log_P(f.ab, idx), 'ab');
subplot(5,2,8); plot_one(t(idx), log_X(f.g0, idx), log_P(f.g0, idx), 'g0');
subplot(5,2,9); plot_one(t(idx), log_X(f.rbc, idx), log_P(f.rbc, idx), 'rbc');
subplot(5,2,10); plot_one(t(idx), log_X(f.tbc, idx), log_P(f.tbc, idx), 'tbc');
end

function [] = plot_one(t, x, p, name)
plot(t, x(1,:), 'r',   t, x(2,:), 'g',   t, x(2,:), 'b',...
     t, p(1,:), 'r-.', t, p(2,:), 'g-.', t, p(3,:), 'b-.',...
     t,-p(1,:), 'r-.', t,-p(2,:), 'g-.', t,-p(3,:), 'b-.'); title(name);
 
end

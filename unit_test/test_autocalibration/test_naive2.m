% Naive test2: Test IMU Kalman Filter
% Tasks:
% 1. We know the starting states (except angular acceleration), given the
% imu measurements, can we reconstruct the trajectory?
% 2. How does the covariance of the states evolve?
% 2. Can I add gravity into the state?
% 3. Can I add sensor bias into the state?

function [f] = test_naive2()

% Add path
addpath('../../lie_group');

% Load simulation environment
load('env_naive.mat');

%% To get comparable at, we need to process it first:
t = full_states.t;
at_gt = zeros(3, numel(t));
for i=1:numel(t);
   R = screw_exp(full_states.rx(:,i));
   at_gt(:,i) = R*full_states.at(:,i);
end

%% Initialize the filter
X_init = zeros(24,1);
X_init( 1: 3) = full_states.rx(:,1);
X_init( 4: 6) = full_states.tx(:,1);
% X_init( 7: 9) = full_states.wt(:,1);
R_im = screw_exp(X_init(1:3));
X_init(10:12) = R_im*full_states.vt(:,1);
%X_init(16:18) = R_im*full_states.at(:,1); % TODO:
X_init(19:21) = [pi 0 0]';
X_init(22:24) = [0 0 0]';

P_rx = 0.0*eye(3);
P_cx = 0.0*eye(3);
P_wt = 10.0*eye(3);
P_vt = 0.0*eye(3);
P_at = 10.0*eye(3);
P_nt = 0.01^2*eye(3);
P_rxic = 0.02^2*eye(3);
P_cxic = 0.02^2*eye(3);

P_init = blkdiag(P_rx, P_cx, P_wt, P_vt, P_nt, P_at, P_rxic, P_cxic);

f = autocalibration_naive(X_init, P_init);
clear -regexp ^P_; clear X_init;

%% Initialize logging variables
num_meas = numel(measurement.TimeStamp);
imu_meas_idx = measurement.Type==0;
t_meas = measurement.TimeStamp(imu_meas_idx);
num_imu_meas = numel(t_meas);
rx_it = zeros(3, num_imu_meas); rx_p = zeros(3, num_imu_meas);
tx_it = zeros(3, num_imu_meas); tx_p = zeros(3, num_imu_meas);
wt_it = zeros(3, num_imu_meas); wt_p = zeros(3, num_imu_meas);
at_it = zeros(3, num_imu_meas); at_p = zeros(3, num_imu_meas);
j = 1;
%% Start the simulation
for i=1:num_meas
    if measurement.Type(i) == 0 % IMU update:
        w = measurement.Data(i, 1:3)';
        a = measurement.Data(i, 4:6)';
        R_imu = 1e-8*eye(6);

        f = f.onImuUpdate(w, a, R_imu, measurement.TimeStamp(i));
        
        p = diag(f.P);
        rx_it(:, j) = f.X(f.rx); rx_p(:, j) = p(f.rx);
        tx_it(:, j) = f.X(f.tx); tx_p(:, j) = p(f.tx);
        wt_it(:, j) = f.X(f.wt); wt_p(:, j) = p(f.wt);
        at_it(:, j) = f.X(f.at); at_p(:, j) = p(f.at);
        
        j = j+1;
    end
end

figure(2); clf;
subplot(2,2,1); plot_data(t, full_states.rx, t_meas, rx_it); title('rx');
subplot(2,2,2); plot1(t_meas, rx_p);

subplot(2,2,3); plot_data(t, full_states.tx, t_meas, tx_it); title('tx');
subplot(2,2,4); plot1(t_meas, tx_p);

% figure(2); clf;
% subplot(2,2,1); plot_data(t, full_states.rx, t_meas, rx_it); title('rx');
% subplot(2,2,2); plot_data(t, full_states.tx, t_meas, tx_it); title('tx');

% subplot(2,2,3); plot_data(t, full_states.wt, t_meas, wt_it); title('wt');
% subplot(2,2,4); plot_data(t, at_gt, t_meas, at_it);          title('at');
end

function [] = plot_data(t1, x1, t2, x2)
plot(t1, x1(1,:), 'r-', t1, x1(2,:), 'g-', t1, x1(3,:), 'b-'); hold on;
plot(t2, x2(1,:), 'r-.', t2, x2(2,:), 'g-.', t2, x2(3,:), 'b-.'); hold off;
end

function [] = plot1(t1, x1)
plot(t1, x1(1,:), 'r-', t1, x1(2,:), 'g-', t1, x1(3,:), 'b-');
end
% Naive test1: Numerical integration
% Tasks:

function test_naive
% Add path
addpath('../../lie_group');

% Load simulation environment
load('env_naive.mat');

% % Initialize the filter
% X_init = zeros(24, 1);
% P_init = 100*eye(24);
% filter = autocalibration_naive(X, P);

% For all time steps, run the filter
t = full_states.t;
N = numel(t);

dt = t(2) - t(1);

rx_it = zeros(3, N);
tx_it = zeros(3, N);
wt_it = zeros(3, N);
vt_it = zeros(3, N);
at_it = zeros(3, N); at_gt = zeros(3, N);

rx_it(:,1) = full_states.rx(:,1);
tx_it(:,1) = full_states.tx(:,1);
R_im = screw_exp(rx_it(:,1));
vt_it(:,1) = R_im*full_states.vt(:,1);
at_it(:,1) = R_im*full_states.at(:,1);
at_gt(:,1) = R_im*full_states.at(:,1);
wt_it(:,1) =      full_states.wt(:,1);
for i = 2:N
    
    R_im = R_im*screw_exp(wt_it(:,i-1)*dt);
    rx_it(:,i) = screw_log(R_im);
    tx_it(:,i) = tx_it(:,i-1) + vt_it(:,i-1)*dt;
    
    vt_it(:,i) = vt_it(:,i-1) + at_it(:,i-1)*dt;
    at_it(:,i) = R_im * full_states.at(:,i);
    at_gt(:,i) = screw_exp(full_states.rx(:,i))*full_states.at(:,i);
    wt_it(:,i) = full_states.wt(:,i);
end

%vt_gt = [-sin(full_states.a); zeros(size(t)); cos(full_states.at)]*full_states.ad;

figure(1); clf;
subplot(2,2,1); plot_data(t, full_states.rx, t, rx_it); title('rx');
subplot(2,2,2); plot_data(t, full_states.tx, t, tx_it); title('tx');

subplot(2,2,3); plot_data(t, full_states.wt, t, wt_it); title('wt');
subplot(2,2,4); plot_data(t, at_gt, t, at_it);          title('at');
end

function [] = plot_data(t1, x1, t2, x2)
plot(t1, x1(1,:), 'r-', t1, x1(2,:), 'g-', t1, x1(3,:), 'b-'); hold on;
plot(t2, x2(1,:), 'r-.', t2, x2(2,:), 'g-.', t2, x2(3,:), 'b-.'); hold off;
end
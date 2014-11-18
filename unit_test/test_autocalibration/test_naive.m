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

rx_it(:,1) = full_states.rx(:,1);
tx_it(:,1) = full_states.tx(:,1);
R_im = screw_exp(rx_it(:,1));
vt_it = R_im*full_states.vt(:,1);
at_it = R_im*full_states.at(:,1);
wt_it =      full_states.wt(:,1);
for i = 2:N
    
    R_im = R_im*screw_exp(wt_it*dt);
    rx_it(:,i) = screw_log(R_im);
    tx_it(:,i) = tx_it(:,i-1) + vt_it*dt;
    
    vt_it = vt_it + at_it * dt;
    at_it = R_im * full_states.at(:,i);
    wt_it = full_states.wt(:,i);
end

%vt_gt = [-sin(full_states.a); zeros(size(t)); cos(full_states.at)]*full_states.ad;

subplot(2,1,1);
plot(t, full_states.rx(1,:), 'r', t, full_states.rx(2,:), 'g', t, full_states.rx(3,:), 'b');
hold on;
plot(t, rx_it(1,:), 'r-.', t, rx_it(2,:), 'g-.', t, rx_it(3,:), 'b-.');
hold off;
subplot(2,1,2);
plot(t, full_states.tx(1,:), 'r', t, full_states.tx(2,:), 'g', t, full_states.tx(3,:), 'b');
hold on;
plot(t, tx_it(1,:), 'r-.', t, tx_it(2,:), 'g-.', t, tx_it(3,:), 'b-.');
hold off;
% subplot(3,1,3);
% plot(t, vt_gt(1,:), 'r', t, vt_gt(2,:), 'g', t, vt_gt(3,:), 'b');
% hold on;
% plot(t, tx_it(1,:), 'r-.', t, tx_it(2,:), 'g-.', t, tx_it(3,:), 'b-.');
% hold off;
end
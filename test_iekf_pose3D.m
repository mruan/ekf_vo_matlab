%clear;
rng(0);

addpath('lie_group');

% Load simulation data:
load('simulation/sim_circle.mat');

j = 100;

sdv_r = 0.05;
sdv_t = 0.05;
sigma = blkdiag(sdv_r^2*eye(3), sdv_t^2*eye(3));
delta_x = sqrt(sigma)*randn(6,1);

gt_rx = trajectory.Rx(:,j);
gt_cx = trajectory.Tx(:,j); 

rx0 = screw_log(screw_exp(delta_x(1:3))*screw_exp(gt_rx));
cx0 = delta_x(4:6) + gt_cx; 

x3d = features.x3d(:, observes.tag{j});
NumPts = size(x3d, 2);
% assume we can measure the 3D points directly:
x2d = observes.x2d{j} + 1*randn(2, NumPts);

filter = iekf_pose_2d([gt_rx; gt_cx], sigma);

filter = filter.update(camera, x2d, x3d);

[gt_rx' gt_cx'; rx0' cx0'; filter.X']
sqrt(abs(filter.P))

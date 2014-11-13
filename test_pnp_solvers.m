% Test if iterative Lie PnP solver actually works

%clear;
rng(1);

addpath('lie_group');
addpath('epnp');

% Load simulation data:
load('simulation/sim_01_cardioid/sim_01.mat');

j = 1;

sdv_r = 0.05;
sdv_t = 0.05;
sigma = blkdiag(sdv_r^2*eye(3), sdv_t^2*eye(3));
delta_x = sqrt(sigma)*randn(6,1);

gt_rx = screw_log(trajectory.R(:,:,j));
gt_cx = trajectory.t(:, j); 

rx0 = screw_log(screw_exp(delta_x(1:3))*trajectory.R(:,:,j));
cx0 = delta_x(4:6) + gt_cx; 

x3d = features.where(:, measurements.feature_tags{j});
NumPts = size(x3d, 2);
% assume we can measure the 3D points directly:
x2d = measurements.image_coords{j} + 0.1*randn(2, NumPts);

x3d_h = [x3d; ones(1, NumPts)]';
x2d_h = [x2d; ones(1, NumPts)]';


% Solve PnP with EPnP, EPnP-GN and PnP-Lie-GN
[R_epnp, T_epnp, ~, ~] = efficient_pnp(x3d_h, x2d_h, camera.K);
rx_epnp = screw_log(R_epnp);

[R_egn, T_egn, ~, ~] = efficient_pnp_gauss(x3d_h, x2d_h, camera.K);
rx_egn = screw_log(R_egn);

[rx_lie, T_lie, flag] = pnp_LieGN(x3d, x2d, camera.K, rx0, cx0);

% Show and compare the results:
[gt_rx'; rx0'; rx_epnp'; rx_egn'; rx_lie']
[gt_cx'; cx0'; T_epnp'; T_egn'; T_lie']


% Conclusion:
% 1. All converge to ground truth without noise
% 2. With noise, epnp tends to drift away. lieGN is better than epnp-gn
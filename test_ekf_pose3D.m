%clear;
rng(0);

addpath('lie_group');

% Load simulation data:
load('simulation/sim_01_cardioid/sim_01.mat');

j = 1;

sdv_r = 0.05;
sdv_t = 0.05;
sigma = blkdiag(sdv_r^2*eye(3), sdv_t^2*eye(3));
delta_x = sqrt(sigma)*randn(6,1);

rx = screw_log(screw_exp(delta_x(1:3))*trajectory.R(:,:,j));
cx = delta_x(4:6) + trajectory.t(:, j); 

filter = ekf_pose_3D([rx; cx], sigma);
ftr_j = features.where(:, measurements.feature_tags{j});
% assume we can measure the 3D points directly:
obs_j = features.where(:, measurements.feature_tags{j});
% then add some Gaussian noise
obs_j = 0.02*randn(3, size(obs_j,2)) + obs_j;

f = filter.update(obs_j(:,1), ftr_j(:,1));

% Print output:
filter.P
f.P

filter.X
f.X


%{
gt = twist_log([trajectory.R(:,:,j) trajectory.t(:, j); 0 0 0 1])'

gt = gt'+ delta_x;
%}
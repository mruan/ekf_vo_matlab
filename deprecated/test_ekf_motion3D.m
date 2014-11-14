%clear;
rng(0);

addpath('lie_group');

% Load simulation data:
load('simulation/sim_01_cardioid/sim_01.mat');

j = 1;

sigma = blkdiag(0.05*eye(3), 0.05*eye(3));
delta_x = sigma*randn(6,1);

rx = screw_log(screw_exp(delta_x(1:3))*trajectory.R(:,:,j));
cx = delta_x(4:6) + trajectory.t(:, j); 

gt = twist_log([trajectory.R(:,:,j) trajectory.t(:, j); 0 0 0 1])'

gt = gt'+ delta_x;
filter = ekf_motion_3D([rx; cx], eye(6));

filter.X'
%filter.P

obs_j = measurements.image_coords{j};
ftr_j = features.where(:, measurements.feature_tags{j});
%decoupled_
filter = filter.update(camera, obs_j(:,:), ftr_j(:,:));

filter.X'
%filter.P

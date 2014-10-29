clear;
rng(0);

addpath('lie_group');

% Load simulation data:
load('simulation/sim_01_cardioid/sim_01.mat');

j = 1;

sigma = 0.01*eye(6);
delta_x = sigma*randn(6,1);

rx = screw_log(screw_exp(delta_x(1:3))*trajectory.R(:,:,j));
cx = delta_x(4:6) + trajectory.t(:, j);

filter = ekf_motion_3D([rx; cx], sigma);

filter.X
filter.P

obs_j = measurements.image_coords{j};
ftr_j = features.where(:, measurements.feature_tags{j});

filter.update(camera, obs_j, ftr_j);

filter.X
filter.P



% Repeatable randomness
rng(1);

% Setup the camera:
camera = generate_camera_info();

% Populate the environment with features:
% z-x defines the plane, y defines the thickness
bb = [-60 60; -15 15; -60 40];
num_features = 100;
features = generate_random_features(bb, num_features);

% Generate camera trajectory:
radius = 10;
trajectory = generate_cardoid_trajectory(radius);

%%
plot_sim_env(trajectory, features);
% Get the visible set of features for each frame:
measurements = generate_measurements(camera, features, trajectory);

save('sim_01.mat', 'camera', 'features', 'measurements', 'trajectory');
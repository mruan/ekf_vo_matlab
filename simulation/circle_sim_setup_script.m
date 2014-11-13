
[camera, trajectory, features, observes] = circle_sim_environment();

plot_sim_env(trajectory, features);

save('sim_circle.mat', 'camera', 'features', 'observes', 'trajectory');
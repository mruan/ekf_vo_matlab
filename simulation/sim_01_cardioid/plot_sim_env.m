function [] = plot_sim_env( trajectory, features )
%PLOT_SIM_ENV Summary of this function goes here
%   Detailed explanation goes here

figure(1); clf; hold on;

N_timestamps = size(trajectory.t, 2);
%N_features   = size(features.tag, 2);

for i= 1:20:N_timestamps
    plot_traj(trajectory.R(:,:,i), trajectory.t(:,i));
end

plot3(features.where(1,:), features.where(2,:), features.where(3,:), '*b');
hold off;

xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; 
end

function [] = plot_traj(R, t)

C = -R'*t;

sc = 0.75;
draw_axis(C, C+sc*R(1,:)', 'r');
draw_axis(C, C+sc*R(2,:)', 'g');
draw_axis(C, C+sc*R(3,:)', 'b');

end

function [] = draw_axis(s, e, spec)

plot3([s(1) e(1)], [s(2) e(2)], [s(3) e(3)], spec);
end
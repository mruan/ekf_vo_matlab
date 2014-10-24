% clear all;

% Where is the target?
target = [2 13 0]';
look_at= [0 15 0]';
up_dir = [0 0 -1]';
radius = 15;

N = 10;
M = traj_3D_setup(target, look_at, up_dir, radius, N);

width = 640;
height= 480;
K = [300 0 width/2; 0 300 height/2; 0 0 1]; % Camera intrinsic matrix

% Initialize the Kalman Filter
filter = EKF_3D([2.5, 13.5, 0.5]', diag([10, 10, 10]));

figure(2); clf; axis([0, width, 0 height]); hold on;
for i=1:N
    P = M(:,1:3,i)*target+M(:,4,i);
    p = K*P;

    x = [p(1)/p(3), p(2)/p(3)]';
%%  Debug code: Show the projected point on the image
    u = p(1)/p(3);
    v = height - p(2)/p(3);
    scatter(u, v, '*r');
    clear u v;  % clear variables afterwards

%%  Add Gaussian Noise to it:


%%  Run the filter
    Q_noise = diag([0.00, 0.00, 0.00]);    % Unit in world scale
    R_noise = diag([0.5, 0.5]);         % Unit in pixel
    filter = filter.update(x, K*M(:,:,i), Q_noise, R_noise);
    filter.X
    filter.P
    
%%  Plot updated projection
    P = M(:,1:3,i)*filter.X+M(:,4,i);
    p = K*P;
    x = [p(1)/p(3), p(2)/p(3)]';
    u = p(1)/p(3);
    v = height - p(2)/p(3);
    scatter(u, v, '*b');
    clear u v;  % clear variables afterwards
end 
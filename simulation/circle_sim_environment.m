function [cam, traj, features, observes] = circle_sim_environment()
%CIRCLE_SIM_ENVIRONMENT Summary of this function goes here
%   Detailed explanation goes here

% Parameters that control the trajectory
% A circle in the xz-plane, centered at (radius, 0 0)
traj_radius  = 10;
NumTrajSteps = 1000;

% Parameters that control the feature placement
feature_radius = 15; % the features are on a bigger circle
NumCols = 50;        % where to place the features
ColHeight = 5;       % range of where to put the features
NumFeatPerCol = 11;  % how many per column

cam  = generate_camera_info();

traj = generate_circle_trajectory();

features = generate_features();

observes = generate_observes(); % use previous outputs: cam, traj, features

%% End of body

function [camera] = generate_camera_info()

% A simple pin-hole camera model

camera.K = [500 0 320; 0 500 240; 0 0 1];
camera.width = 640;
camera.height= 480;

end

%% Generate circle trajectory
function [trajectory] = generate_circle_trajectory()
    theta = linspace(0, 2*pi, NumTrajSteps);

    Tx = zeros(3, NumTrajSteps);
    Rx = zeros(3, NumTrajSteps);

    for i=1:NumTrajSteps
        a = theta(i);

        Tx(:, i) = [traj_radius*(1-cos(a)); 0.0; traj_radius*sin(a)];

        R = [cos(a) 0 sin(a);...
            0  1     0 ;...
            -sin(a) 0 cos(a)];
        Rx(:, i) = screw_log(R);
    end

    trajectory.Tx = Tx;
    trajectory.Rx = Rx;
end

%% Generate features
function [features] = generate_features()

% coincide with traj center
center_x = traj_radius;
center_z = 0.0;

theta = linspace(0, 2*pi, NumCols);

% Pre-allocate space
features.x3d = zeros(3, NumCols*NumFeatPerCol);
features.tag = 1:NumCols*NumFeatPerCol;

y_row = linspace(-ColHeight, ColHeight, NumFeatPerCol);
base = 1;
for i=1:NumCols
    a = theta(i);
    x_row = repmat(center_x - feature_radius*cos(a), 1, NumFeatPerCol);
    z_row = repmat(center_z + feature_radius*sin(a), 1, NumFeatPerCol);
    
    features.x3d(:,base:base+NumFeatPerCol-1) = [x_row; y_row; z_row];
    base = base + NumFeatPerCol;
end

end

%% Generate observations
function [observes] = generate_observes()

K = cam.K;
W = cam.width;
H = cam.height;

N_poses    = size(traj.Tx, 2);
N_features = size(features.x3d, 2);

observes.x2d = cell(N_poses, 1);
observes.tag = cell(N_poses, 1);
for i=1:N_poses
   % Convert [Rs Ts] (camera-in-world) -> [Rc Tc] (world-in-camera)
   Rc = screw_exp(traj.Rx(:, i))';
   Tc = -Rc*traj.Tx(:, i);
   
   cam_coord = Rc*features.x3d + repmat(Tc, 1, N_features);
   
   % Try to project each feature into the current frame:
   % 1. Are they in front of the camera?
   not_in_front = cam_coord(3, :) < 0;
   
   img_coord = K*cam_coord;
   img_coord = [img_coord(1,:)./img_coord(3,:); img_coord(2,:)./img_coord(3,:)];
   
   % 2. Are they visible in the frame?
   not_in_frame = img_coord(1,:) < 0 | img_coord(1,:) > W | img_coord(2,:) < 0 | img_coord(2,:) > H;
   
   % 3. Combine 1 and 2:
   not_visible = not_in_front | not_in_frame;
   
   observes.x2d{i} = img_coord(:, ~not_visible);
   observes.tag{i} = features.tag(~not_visible);
   
   %% Verbose: print out some debug warnings:
   num_pts_in_view = numel(observes.x2d{i});
   min_pts_in_view = 5;
   if (num_pts_in_view < min_pts_in_view)
       fprintf('Frame %d has  %d features (min of %d required)\n', i, num_pts_in_view, min_pts_in_view);
   end
end
end

end



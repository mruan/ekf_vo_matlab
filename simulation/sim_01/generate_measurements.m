function [measurement] = generate_measurements(camera, features, trajectory)
%MEASUREMENT_SETUP Summary of this function goes here
%   Detailed explanation goes here

K = camera.K;
W = camera.width;
H = camera.height;

N_timestamps = size(trajectory.t, 2);
N_features   = size(features.tag, 2);

measurement.image_coords = cell(N_timestamps, 1); % where the feature is projected onto the image?
measurement.feature_tags = cell(N_timestamps, 1); % which feature does the imaged point belong to?

% For each time stamp
for i=1:N_timestamps
    
    % Try to project each feature into current frame:
    cam_coord = trajectory.R(:,:,i)*features.where+repmat(trajectory.t(:,i), 1,N_features);

    % Is the feature in front of the camera?
    not_in_front = cam_coord(3,:) < 0;
    
    img_coord = K*cam_coord;
    img_coord = [img_coord(1,:)./img_coord(3,:); img_coord(2,:)./img_coord(3,:)];

    % Is the feature visible in the image?
    not_in_frame = img_coord(1,:) < 0 | img_coord(1,:) > W | img_coord(2,:) < 0 | img_coord(2,:) > H;
    
    % Not visible if either of the checks fails
    not_visible = not_in_front | not_in_frame;
    
    measurement.feature_tags{i} = features.tag(~not_visible);
    measurement.image_coords{i} = img_coord(:, ~not_visible);
    
    %% Verbose: print out some debug info:
    num_features = numel(measurement.feature_tags{i});
    min_features = 5;
    if (num_features < min_features)
       fprintf('Frame %d has less than %d features in the frame: %d\n', i, min_features,  num_features);
    end
end

end


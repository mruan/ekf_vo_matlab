function [features] = generate_random_features(bb, Num_features)

% Populate the environment with features

% min_x = bb(1,1);
% max_x = bb(1,2);
% min_y = bb(2,1);
% max_y = bb(2,2);
% min_z = bb(3,1);
% max_z = bb(3,2);

% features.where = [(max_x-min_x)*rand(1, Num_features)+min_x;...
%                   (max_y-min_y)*rand(1, Num_features)+min_y;...
%                   (max_z-min_z)*rand(1, Num_features)+min_z];

      
features.tag = 1:Num_features;
radius = 0.5*max(range(bb,2));


min_y = bb(2,1);
max_y = bb(2,2);
% Generate features so that when projected to the z-x plane will look like
% a circle:
% z = radius*cos(t) - offset;
% x = radius*sin(t);
% y = (max_y - min_y)*rand(1, Num_features) + min_y;

t = linspace(0, 2*pi, Num_features);
features.where = [radius*sin(t); ...
                  (max_y - min_y)*rand(1, Num_features) + min_y; ...
                  radius*cos(t)];

end
function [ traj ] = generate_cardoid_trajectory(radius )
%TRAJECTORY_SETUP Summary of this function goes here
%   Detailed explanation goes here

N = 1000;

t = linspace(0, 2*pi, N);

X = radius*(2*sin(t) - sin(2*t));
Y = zeros(size(t));
Z = radius*(2*cos(t) - cos(2*t)-1);

% Orientation 
% special cases of 1 and N where the time derivative is 0;
R(:,:,1) = [1 0 0; 0 1 0; 0 0 1];
R(:,:,N) = [-1 0 0; 0 -1 0; 0 0 -1];

% the T (translation) part
T = zeros(3,N);

for i=2:N-1
    ti = t(i);
    z = [cos(ti)-cos(2*ti); 0; -sin(ti)+sin(2*ti)]/sqrt(2-2*cos(t(i)));
    y = [0 1 0]';
    x = cross(y,z);
   R(:,:,i) = [x y z]';
   
   T(:,i) = -R(:,:,i)*[X(i); Y(i); Z(i)]; % so that C = -R'*t
end

% Pack them up
traj.R = R;
traj.t = T;

end
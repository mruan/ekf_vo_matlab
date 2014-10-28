function [camera] = generate_camera_info()

% A simple pin-hole camera model

camera.K = [500 0 320; 0 500 240; 0 0 1];
camera.width = 640;
camera.height= 480;

end
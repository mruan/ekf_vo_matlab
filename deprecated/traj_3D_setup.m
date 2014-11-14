function [M] = traj_3D_setup(target, look_at, up_dir_, radius_, N)
% Generate a fancy camera trajectory

C = zeros(3, N);
R = zeros(3, 3, N);
M = zeros(3, 4, N);
for i=1:N
   theta = pi/2/N*(i-1);
   
   x = radius_*sin(theta);
   y =  radius_*(1-cos(theta));
   z = 0;
   
   C(:,i)  = [x, y, z]';
   image_y = up_dir_;
   image_z = look_at - [x, y, z]';
   image_z = image_z / norm(image_z);
   image_x = cross(image_y, image_z);
   
   Rc = [image_x image_y image_z];
   
   R(:,:,i) = Rc;
   C(:,i)   = [x, y, z]';
   M(:,:,i) = [Rc' -Rc'*C(:,i)];
end

figure(1); clf; hold on;
plot3(C(1,:), C(2,:), C(3,:), '.r');
scatter3(target(1), target(2), target(3));
for i=1:N
   c  = C(:,i);
   r  = R(:,:,i);
   sc = 0.25;   % control the length of the axes to be displayed

   plot3([c(1) c(1)+sc*r(1,1)], [c(2) c(2)+sc*r(2,1)], [c(3) c(3)+sc*r(3,1)], 'r');
   plot3([c(1) c(1)+sc*r(1,2)], [c(2) c(2)+sc*r(2,2)], [c(3) c(3)+sc*r(3,2)], 'g');
   plot3([c(1) c(1)+sc*r(1,3)], [c(2) c(2)+sc*r(2,3)], [c(3) c(3)+sc*r(3,3)], 'b');
end
hold off;
axis equal;
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

end
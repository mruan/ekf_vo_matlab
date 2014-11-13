function [rx, cx, flag ] = pnp_LieGN( x3d_h, x2d_h, camera_K, rx0, cx0)
%PNP_LIEGN Summary of this function goes here
%   Detailed explanation goes here

max_iter = 200;
min_norm = 1e-7;
converge_flag = false;

N = size(x3d_h, 2);

z = x2d_h(1:2, :);
z = z(:);

rx = rx0; cx = cx0;
for i=1:max_iter
    
    h = zeros(2*N, 1);
    H = zeros(2*N, 6);
    
    fR = screw_exp(rx);
    ft = cx;
    for j = 1:N
        xyz = x3d_h(1:3, j);
        
        uvw = camera_K*(fR*xyz + ft);
        
        uv  = uvw(1:2)./uvw(3);
        
        h(2*j-1 :2*j) = uv;
        
        H1 = [1/uvw(3)   0.0    -uv(1)/uvw(3);...
            0.0    1/uvw(3) -uv(2)/uvw(3)];
        
        H2 = camera_K*[so3_alg(-fR*xyz) eye(3)];
        
        H(2*j-1 :2*j, :) = H1*H2;
    end
    
    y = z - h;
    dx = (H'*H)\H'*y;
    
    if norm(dx) < min_norm
        converge_flag = 1;
        break;
    end
    
    rx = screw_log(screw_exp(dx(1:3))*screw_exp(rx));
    cx = cx + dx(4:6);
end

if converge_flag
    fprintf('Converged after %d iterations\n', i);
else
    fprintf('Failed to converged after %d iterations\n', max_iter);
end

flag = converge_flag;
end



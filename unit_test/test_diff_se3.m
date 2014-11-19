function [] = test_diff_se3()
rng(1);

K = [500, 0, 320; 0, 500, 240; 0, 0, 1];

xyz = 20*rand(3,1);

rx = rand(3,1);
cx = rand(3,1);

Rx = screw_exp(rx);

% uvw = img_proj_bs(Rx, cx, xyz);
% uv = img_proj_bs_norm(Rx, cx, xyz);

uvw = img_proj_sb(Rx, cx, xyz);
uv =  img_proj_sb_norm(Rx, cx, xyz);


delta = 1e-6;
Hn = zeros(2,6);
for i=1:3
    d_rx = [0 0 0]';
    d_rx(i) = delta;
    
    Rx_p = screw_exp(d_rx)*Rx;
    
%     uvw_p = img_proj(Rx_p, cx, xyz);
%     Hn(:,i) = (uvw_p-uvw)/delta;

%     uv_p = img_proj_bs_norm(Rx_p, cx, xyz);
    uv_p = img_proj_sb_norm(Rx_p, cx, xyz);
    
    Hn(:,i) = (uv_p-uv)/delta;
end

for i=4:6
    d_rx = [0 0 0]';
    d_rx(i-3) = delta;
    
%     uvw_p = img_proj(Rx, cx+d_rx, xyz);
%     Hn(:,i) = (uvw_p-uvw)/delta;
    
%     uv_p = img_proj_bs_norm(Rx, cx+d_rx, xyz);
    uv_p = img_proj_sb_norm(Rx, cx+d_rx, xyz);
    
    Hn(:,i) = (uv_p-uv)/delta;
end
% H numerical
Hn 

% H analytical
% Ha = [1/uvw(3) 0 -uv(1)/uvw(3);...
%       0 1/uvw(3) -uv(2)/uvw(3)]*K*[so3_alg(-Rx*xyz) eye(3)]
Ha = [1/uvw(3) 0 -uv(1)/uvw(3);...
      0 1/uvw(3) -uv(2)/uvw(3)]*K*Rx'*[so3_alg(xyz-cx) -eye(3)]
  
    function [uvw] = img_proj_bs(R, t, p)       
        uvw = K*(R*p+t);
    end

    function [uvw] = img_proj_sb(R, t, p)
        uvw = K*R'*(p-t);
    end

    function [uv] = img_proj_bs_norm(R, t, p)
       uv_un = img_proj_bs(R, t, p);
       uv = uv_un(1:2)./uv_un(3);
    end

    function [uv] = img_proj_sb_norm(R, t, p)
       uv_un = img_proj_sb(R, t, p);
       uv = uv_un(1:2)./uv_un(3);
    end
end


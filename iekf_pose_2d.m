classdef iekf_pose_2d
    %IEKF_POSE_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X; % [rx cx] as body-in-world
        P;
        
        max_iteration;
        min_variation;
    end
    
    methods
        function f = iekf_pose_2d(varargin)
            if nargin == 0
                f.X = zeros(6, 1);
                %f.P = blkdiag(zeros(6), 500*ones(6));
                f.P = eye(6);
            elseif nargin == 2
                f.X = varargin{1};
                f.P = varargin{2};
            end   
            f.max_iteration = 50;
            f.min_variation = 1e-5;
        end
        
        function f = update(f, camera, observes, features)
           num_pts = size(observes, 2);

           z = observes(:);
           R = 1.0*eye(2*num_pts);  % image measurement noise
           
           Rx0 = screw_exp(f.X(1:3));
           Tx0 = f.X(4:6);
           P0  = f.P;

           Rx = Rx0;
           Tx = Tx0;
           dx_im1 = zeros(6,1);
           for i = 1:f.max_iteration
               
               [hi, Hi] = inner_loop();
               
               yi = z - hi + Hi*dx_im1; % IEKF innovation
               Si = Hi*P0*Hi' + R;    
               Ki = P0*Hi'/Si;
               
               dx_i = Ki*yi;
               
               f.X(1:3) = screw_log(screw_exp(dx_i(1:3))*Rx0);
               f.X(4:6) = dx_i(4:6) + Tx0;
           
               if (norm(dx_i - dx_im1) < f.min_variation)
                   f.P = P0 - Ki*Hi*P0;
                   break;
               else
                   dx_im1 = dx_i;
               end
               
               % Update current estimate of the states
               Rx = screw_exp(f.X(1:3));
               Tx = f.X(4:6);
           end
           
            function [h, H] = inner_loop()
                h = zeros(2*num_pts, 1);
                H = zeros(2*num_pts, 6);
                
                for j=1:num_pts
                    x3d_s = features(:,j);    % pts in inertial frame
                    x3d_b = Rx'*(x3d_s - Tx); % pts in body frame
                    
                    x2d_h = camera.K * x3d_b; % pts in image homo. frame
                    x2d_i = [x2d_h(1)./x2d_h(3); x2d_h(2)./x2d_h(3)];
                    
                    H1 = [1/x2d_h(3)   0         -x2d_i(1)/x2d_h(3);...
                            0        1/x2d_h(3)  -x2d_i(2)/x2d_h(3)];
                    
                    H2 = camera.K*Rx'*[so3_alg(x3d_s-Tx) -eye(3)];
                    
                    h(2*j-1:2*j)   = x2d_i;
                    H(2*j-1:2*j,:) = H1*H2;
                end
            end
        end
        
    end
    
end


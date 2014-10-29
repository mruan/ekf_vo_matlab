classdef ekf_motion_3D
    %EKF_MOTION_3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X; %[rot(x y z) (x y z)] % ang(vx vy vz) (vx vy vz)]
        P; % 6x6                 % 12x12 covariance matrix
    end
    
    methods
        function f = ekf_motion_3D(varargin)
            if nargin == 0
                f.X = zeros(6, 1);
                %f.P = blkdiag(zeros(6), 500*ones(6));
                f.P = eye(6);
            elseif nargin == 2
                f.X = varargin{1};
                f.P = varargin{2};
            end
        end
        
  %{      
        function f = predict(f, delta_t, LIN_ACC, ANG_ACC)
            % Constant acceleration model: dx is thus a Gaussian
            mean_delta_vt = f.X(7:12)*delta_t;
            
            % v~N(v0, Sigma_v), a~N(0, Sigma_a)
            % dx = v*dt + 0.5*dt^2*a=[dt, 0.5*dt^2]*[v a]'
        end
  %}      
        function f = update(f, camera, measurements, features)
            
            num_features = size(features, 2);
            
            fR = screw_exp(f.X(1:3));
            ft = f.X(4:6);
            
            z = measurements(:); % rasterize to a tall column vector
            h = zeros(2*num_features, 1);
            H_x = zeros(2*num_features, 3);
            H_r = zeros(2*num_features, 3);
            % TODO: vectorize this loop
            for i=1:num_features
                xyz = features(:, i);
                
                uvw = camera.K*(fR*xyz + ft);
                
                uv  = uvw(1:2)./uvw(3);
                
                h(2*i-1:2*i) = uv;
                
                H1 = [1/uvw(3)   0.0    -uv(1)/uvw(3);...
                        0.0    1/uvw(3) -uv(2)/uvw(3)];
                %  
                H2 = camera.K*[so3_alg(-fR*xyz)];
                
                H_r(2*i-1:2*i,:) = H1*H2;
                H_x(2*i-1:2*i,:) = H1*camera.K;
            end
            
            % From here on R and K refer to the Kalman context
            R = sparse(diag(1.0*ones(2*num_features, 1)));
            
            y = z - h;
            S = H_r*f.P(1:3,1:3)*H_r'+ R;
            K = f.P(1:3,1:3)*H_r'/S;
            dx = K*y;
            
            f.P(1:3,1:3) = (eye(3) - K*H_r)*f.P(1:3,1:3);
            
            % Note the different update rules for the states:
            screw_add =@(dx, x) screw_log(screw_exp(dx)*screw_exp(x));
            f.X(1:3) = screw_add(dx, f.X(1:3));
            
            S = H_x*f.P(4:6,4:6)*H_x'+ R;
            K = f.P(4:6,4:6)*H_x'/S;
            dx = K*y;           
            f.X(4:6) = dx + f.X(4:6);
            
            % Update the covariance matrix
            f.P(4:6,4:6) = (eye(3) - K*H_x)*f.P(4:6,4:6);
        end
    end
    
end


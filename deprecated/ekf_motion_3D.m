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
        function f = decoupled_update(f, camera, meas, feat)
            num_features = size(feat, 2);
            
            fR = screw_exp(f.X(1:3));
            ft = f.X(4:6);
            
            z = meas(:); % rasterize to a tall column vector
            h = zeros(2*num_features, 1);
            Hr = zeros(2*num_features, 6);
            
            % Do rotation first:
            for i=1:num_features
                xyz = feat(:, i);
                
                uvw = camera.K*(fR*xyz + ft);
                
                uv  = uvw(1:2)./uvw(3);
                
                h(2*i-1:2*i) = uv;
                
                H1 = [1/uvw(3)   0.0    -uv(1)/uvw(3);...
                        0.0    1/uvw(3) -uv(2)/uvw(3)];
                
                H2 = camera.K*[so3_alg(-fR*xyz) zeros(3)];
                
                Hr(2*i-1:2*i,:) = H1*H2;
            end
            
            % From here on R and K refer to the Kalman context
            R = sparse(diag(100.0*ones(2*num_features, 1)));
            
            y = z - h;
            S = Hr*f.P*Hr'+ R;
            K = f.P*Hr'/S;
            dx = K*y;
            
            % Note the different update rules for the states:
            f.X(1:3) = screw_log(screw_exp(dx(1:3))*screw_exp(f.X(1:3)));
            Kp = (eye(6) - K*Hr);
            f.P= Kp*f.P*Kp'+K*R*K';
            
%             Hx = zeros(2*num_features, 6);
%             fR = screw_exp(f.X(1:3));
%             
%             for i=1:num_features
%                 xyz = feat(:, i);
%                 
%                 uvw = camera.K*(fR*xyz + ft);
%                 
%                 uv  = uvw(1:2)./uvw(3);
%                 
%                 h(2*i-1:2*i) = uv;
%                 
%                 H1 = [1/uvw(3)   0.0    -uv(1)/uvw(3);...
%                         0.0    1/uvw(3) -uv(2)/uvw(3)];
%                 
%                 H2 = camera.K*[zeros(3) eye(3)];
%                 Hx(2*i-1:2*i,:) = H1*H2;
%             end
%             y = z - h;
%             R = sparse(diag(1.0*ones(2*num_features, 1)));
%             S = Hx*f.P*Hx'+ R;
%             K = f.P*Hx'/S;
%             dx = K*y;
%             
%             f.X(4:6) = dx(4:6) + f.X(4:6);
%             f.P = (eye(6) - K*Hx)*f.P;
        end
        
        function f = update_se3(f, camera, meas, feat)
            num_features = size(feat, 2);
            
            T = twist_exp(f.X);
            K = [camera.K zeros(3,1)];
            
            z = meas(:); % rasterize to a tall column vector
            h = zeros(2*num_features, 1);
            H = zeros(2*num_features, 6);
            
            for i=1:num_features
                xyz1 = [feat(:, i); 1];
                
                t_xyz1 = T*xyz1;
                uvw = K*t_xyz1;
                
                uv  = uvw(1:2)./uvw(3);
                
                h(2*i-1:2*i) = uv;
                
                H1 = [1/uvw(3)   0.0    -uv(1)/uvw(3);...
                        0.0    1/uvw(3) -uv(2)/uvw(3)];
                %so3_alg(-t_xyz1(1:3))
                H2 = [zeros(3) eye(3)];
                
                H(2*i-1:2*i,:) = H1*camera.K*H2;
            end
            
            R = sparse(diag(100.0*ones(2*num_features, 1)));
            
            y = z - h;
            S = H*f.P*H'+ R;
            K = f.P*H'/S;
            dx = K*y;
            
            f.X = twist_log(twist_exp(dx)*twist_exp(f.X));
            Kp = (eye(6) - K*H);
            f.P = Kp*f.P*Kp'+K*R*K';
        end
        
        function f = update(f, camera, meas, feat)
            
            num_features = size(feat, 2);
            
            fR = screw_exp(f.X(1:3));
            ft = f.X(4:6);
            
            z_hat = camera.K\[meas; ones(1,num_features)]; % normalized image coordinate
            z_hat = [z_hat(1,:); z_hat(2,:)];
            z = z_hat(:); % rasterize to a tall column vector
            h = zeros(2*num_features, 1);
            H = zeros(2*num_features, 6);
            % TODO: vectorize this loop
            for i=1:num_features
                xyz = feat(:, i);
                %camera.K*
                uvw = (fR*xyz + ft);
                
                uv  = uvw(1:2)./uvw(3);
                
                h(2*i-1:2*i) = uv;
                
                H1 = [1/uvw(3)   0.0    -uv(1)/uvw(3);...
                        0.0    1/uvw(3) -uv(2)/uvw(3)];
                % camera.K*
                H2 = [so3_alg(-fR*xyz) eye(3)];
                
                H(2*i-1:2*i,:) = H1*H2;
            end
            
            % From here on R and K refer to the Kalman context
            R = sparse(diag(1.0*ones(2*num_features, 1)));
            
            y = z - h;
            S = H*f.P*H'+ R;
            K = f.P*H'/S;
            dx = K*y;
            
            % Note the different update rules for the states:
            screw_add =@(dx, x) screw_log(screw_exp(dx)*screw_exp(x));
            f.X(1:3) = screw_add(dx(1:3), f.X(1:3));
            f.X(4:6) = dx(4:6) + f.X(4:6);
            
            % Update the covariance matrix
            f.P = (eye(6) - K*H)*f.P;
        end
    end
    
end


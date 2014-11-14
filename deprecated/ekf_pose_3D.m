classdef ekf_pose_3D
    %EKF_POSE_3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X; %[rot(x y z) (x y z)
        P; %[6x6]
    end
    
    methods
        
        function f = ekf_pose_3D(varargin)
            if nargin == 0
                f.X = zeros(6, 1);
                %f.P = blkdiag(zeros(6), 500*ones(6));
                f.P = eye(6);
            elseif nargin == 2
                f.X = varargin{1};
                f.P = varargin{2};
            end
        end
        
        function f = update(f, meas, feat)
            num_features = size(feat, 2);
            
            fR = screw_exp(f.X(1:3));
            ft = f.X(4:6);
            
            z = meas(:); % rasterize to a tall column vector
            h = zeros(3*num_features, 1);
            H = zeros(3*num_features, 6);
                       
            for i=1:num_features
                xyz = feat(:, i);
                
                xyz_h = fR*xyz + ft;
                
                h(3*i-2 : 3*i) = xyz_h;
                
                H(3*i-2 : 3*i,:) = [so3_alg(-fR*xyz) eye(3)];
            end
            
            R = sparse(diag(1.0*ones(3*num_features, 1)));
            
            y = z-h;
            S = H*f.P*H';
            S = S + R;
            K = f.P*H'/S;
            dx= K*y;
            
            f.X(1:3) = screw_log(screw_exp(dx(1:3))*screw_exp(f.X(1:3)));
            f.X(4:6) = f.X(4:6) + dx(4:6);
            Kp = (eye(6) - K*H);
            f.P = Kp*f.P*Kp' + K*R*K';
        end
    end
    
end


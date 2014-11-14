classdef ekf_R3
    %EKF_R3 Basically an reimplementation of EKF_3D, except the
    %initialization step is made more general to observe divergence
    
    properties
        x  %[x y z]
        P  % covariance 3x3
        Q  % image space uncertainty: eye(2)
    end
    
    methods
        function obj = ekf_R3()
            obj.x = zeros(3,1);
            obj.P = 1000*eye(3);
            obj.Q = eye(2);
        end
 
        function mu = get.x(obj)
            mu = obj.x;
        end
        
        function cov= get.P(obj)
            cov= obj.P;
        end
        
        % Initialize from depth
        function obj = init(obj, z, K, R, t)
            % Assert uv is a colum vector
            [nr, nc] = size(z);
            assert(nr==2 && nc==1, 'uv must be a col vector');
            
            % unknown depth initialization
            depth_mu    = 10;
            depth_sigma = 3;
            
            F = R'/K;
            obj.x = -R'*t + depth_mu * F *[z; 1];
            A = F *[z; 1];    % A-> 3x1
            B = depth_mu * F;  
            B = B(:,1:2);      % B-> 3x2
            obj.P = A*depth_sigma*A'+B*obj.Q*B';
        end
        
        function obj = update(obj, z, K, R, t)
            % Assert uv is a colum vector
            [nr, nc] = size(z);
            assert(nr==2 && nc==1, 'uv must be a col vector');
            
            % z_h is the homogeneous coord of \hat(z) w/o normalization
            z_h = K*(R*obj.x)+t;
            
            % z_hn i.e. \hat(z) (with normalization in image coord. now)
            z_hn = z_h/z_h(3);
            z_hn = z_hn(1:2);
            
            % Find the linearized observation matrix
            Hp = [1/z_h(3) 0         -z_hn(1)/z_h(3);...
                         0  1/z_h(3) -z_hn(2)/z_h(3)];                   
            H  = Hp*K*R;
            
            %% Do the EKF update
            S  = H*obj.P*H'+obj.Q;
            k  = obj.P*H'/S;
            
            obj.x = obj.x + k*(z - z_hn);
            obj.P = obj.P - k*H*obj.P;
        end
    end
    
end


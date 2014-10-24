classdef ekf_P3
    %EKF_P3 EKF for a point in Projective-3D space
    %   Detailed explanation goes here
    
    properties
        x  %[x y z w]
        P  % covariance 4x4
        Q  % image space uncertainty: eye(2)
    end
    
    methods
        function obj = ekf_P3()
            obj.x = zeros(4,1);
            obj.P = 1000*eye(4);
            obj.Q = eye(2);
        end
        
        function mu = get.x(obj)
            mu = obj.x;
        end
        
        function cov= get.P(obj)
            cov= obj.P;
        end
        
        % Initialize from inverse depth
        function obj = init(obj, z, K, R, t)
            % Assert uv is a colum vector
            [nr, nc] = size(z);
            assert(nr==2 && nc==1, 'uv must be a col vector');
            
            % rho --> inverse depth
            rho_mu    = 1;
            rho_sigma = 1;
            
            F = R'/K;
            obj.x(1:3) = -rho_mu * R'*t + F*[z; 1]; obj.x(4) = rho_mu;
            A = [-R'*t; 1];        % A-> 4x1
            
            B = [F(:,1:2); 0 0];   % B-> 4x2
            obj.P = A*rho_sigma*A'+B*obj.Q*B';
        end
        
        function obj = update(obj, z, K, R, t)
            % Assert uv is a colum vector
            [nr, nc] = size(z);
            assert(nr==2 && nc==1, 'z must be a col vector');
            
            % z_h is the homogeneous coord of \hat(z) w/o normalization
            z_h = K*(R*obj.x(1:3)+t*obj.x(4));
            
            % z_hn i.e. \hat(z) (with normalization in image coord. now)
            z_hn = z_h/z_h(3);
            z_hn = z_hn(1:2);
            
            % Find the linearized observation matrix
            Hp = [1/z_h(3) 0         -z_hn(1)/z_h(3);...
                         0  1/z_h(3) -z_hn(2)/z_h(3)];
            H  = Hp*K*[R t];
            
            %% Do the EKF update
            S  = H*obj.P*H'+obj.Q;
            k  = obj.P*H'/S;
            
            obj.x = obj.x + k*(z - z_hn);
            obj.P = obj.P - k*H*obj.P;
        end
    end
    
end


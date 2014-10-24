classdef EKF_3D
    %EKF_3D A naive impl of EKF to estimate a 3D point
    %       One significant drawback is the requirement for a
    %       decent guess of where the point is (excatly why a 
    %       naive EKF will fail and why EKF_P3 can do better)
    
    properties
        X % Location in 3D
        P % Variance
    end
    
    methods
        function filter_obj = EKF_3D(X_init, P_init)
           % Create a new EKF_3D object
           filter_obj.X = X_init;
           filter_obj.P = P_init;
        end
        
        function mu = get.X(obj)
            mu = obj.X;
        end
        
        function cov= get.P(obj)
            cov= obj.P;
        end
        
        function obj = update(obj, z, M, Q, R)
            P_hat = obj.P + Q;
            
            z_hat = M*[obj.X; 1];         % Predict where X is projected
            z_hat_norm = [z_hat(1)/z_hat(3) z_hat(2)/z_hat(3)]';
            y_hat = z - z_hat_norm;   % compute the residual
            
            
            % H = [fx/Z   0.0  fx*x/z^2] * R
            %     [ 0.0  fy/z  fy*y/z^2]
            H = [1/z_hat(3) 0.0             -z_hat(1)/z_hat(3)^2;
                 0.0             1/z_hat(3) -z_hat(2)/z_hat(3)^2];
            H = H * M(:,1:3);
            S = H*P_hat*H' + R;
            K = P_hat*H'/S;
            
            obj.X = obj.X + K*y_hat;
            obj.P = P_hat - K*H*P_hat;
        end
        
        function plot_err_ellipse(obj)
            % requires the function: error_ellipse from
            % http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse
            % By AJ Johnson
            error_ellipse('mu', obj.X, 'C', obj.P);
        end
    end
    
end


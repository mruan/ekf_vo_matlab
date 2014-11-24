classdef ekf_camcalib
    %EKF_CT_CALIB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    properties(Constant)
       rsb = 1:3;
       tsb = 4:6;
       wt  = 7:9;
       vt  = 10:12;
       at  = 13:15;
       rbc = 16:18;
       tbc = 19:21;
       nonSO3 = [4:15 19:21];
       N_states = 21;
       max_iter = 20;
       dx_threshold = 1e-4;
    end
    
    methods
        function [f] = ekf_camcalib()
%             f.X   = zeros(f.N_states, 1);
%             f.P   = 100*eye(f.N_states);
% 
%             f.time_stamp = 0.0;
        end
        
        %{
        y_imu = [   wt  ]
                [ R'*at ]
        %}
        function [X, P] = onImuUpdate(f, w, a, imu_noise, X_h, P_h)
            % Save variables:
            x  = X_h;
            
            dx = zeros(f.N_states, 1);
            px = dx;
            H  = zeros(6, f.N_states);
            for iter = 1: f.max_iter
               Rsb = screw_exp(x(f.rsb));
               w_pred = x(f.wt);
               a_pred = Rsb'*x(f.at); %
               
               H(1:3, f.wt) = eye(3);
               H(4:6, f.at) = Rsb';
               H(4:6, f.rsb)= Rsb'*so3_alg(x(f.at)); %
              
               y = [w - w_pred; a - a_pred] + H*px;
               
               S = H*P_h*H' + imu_noise;
               K = P_h*H'/S;
               dx= K*y;
               
               x = f.states_add(dx, X_h);
               if max(abs(dx - px)) < f.dx_threshold
                   break;
               end
               px = dx;
            end
            X = x;
            P = P_h - K*H*P_h;
        end

        %{
        y_cam = [  rsc  ]
                [  tsc  ]
        where
        rsc = screw_log(screw_exp(rsb)*screw_exp(rbc))
        tsc = screw_exp(rsb)*tbc + tsb
        %}
        function [X, P] = onCamUpdate(f, rsc, tsc, omni_noise, X_h, P_h)
            % Propagate to current time step
            x = X_h;

            dx = zeros(f.N_states, 1);
            px = dx;
            H  = zeros(6, f.N_states);
            Rsc_obsv = screw_exp(rsc);
            for iter = 1: f.max_iter
               Rsb = screw_exp(x(f.rsb));
               
               Rsc_pred = Rsb*screw_exp(x(f.rbc));
               tsc_pred = Rsb*x(f.tbc) + x(f.tsb);
               
               H(1:3, f.rsb) = eye(3);
               H(1:3, f.rbc) = Rsb;
               H(4:6, f.rsb) = -so3_alg(Rsb*x(f.tbc));
               H(4:6, f.tsb) = eye(3);
               H(4:6, f.tbc) = Rsb;

               y = [screw_log(Rsc_obsv*Rsc_pred'); tsc - tsc_pred] + H*px;
               
               S = H*P_h*H' + omni_noise;
               K = P_h*H'/S;
               dx= K*y;
               
               x = f.states_add(dx, X_h);
               if max(abs(dx - px)) < f.dx_threshold
                   break;
               end
               px = dx;
            end
            X = x;
            P = P_h - K*H*P_h;
        end
        
        function [X_h, P_h] = propagate(f, X, P, dt)
            X_h = X;
            P_h = P;
            % dt could potentially be 0
            if dt < 1e-7
                return;
            end
            
            R = screw_exp(X(f.rsb));
            X_h(f.rsb) = screw_log(R*screw_exp(X(f.wt)*dt));
            X_h(f.tsb) = X(f.tsb) + X(f.vt)*dt;
            X_h(f.vt)  = X(f.vt)  + X(f.at)*dt;
            
            %{
              Linearized model:
              dX = (I+F*dt)*X + (G*dt)*n
            
                    rsb tsb  wt  vt  at rbc tbc
               F =
               rsb [ 0   0   R   0   0   0   0 ]
               tsb [ 0   0   0   I   0   0   0 ]
                wt [ 0   0   0   0   0   0   0 ]
                vt [ 0   0   0   0   I   0   0 ]
                at [ 0   0   0   0   0   0   0 ]
               rbc [ 0   0   0   0   0   0   0 ]
               tbc [ 0   0   0   0   0   0   0 ]
            %}
            F = eye(f.N_states);
            F(f.rsb, f.wt) = R*dt;
            F(f.tsb, f.vt) = eye(3)*dt;
            F(f.vt,  f.at) = eye(3)*dt;
            
            %{         
                    nw  na  nrc ntc
               G =
               rsb [ 0   0   0   0 ]
               tsb [ 0   0   0   0 ]
                wt [ I   0   0   0 ]
                vt [ 0   0   0   0 ]
                at [ 0   I   0   0 ]
               rbc [ 0   0   I   0 ]
               tbc [ 0   0   0   I ]
            %}
            G = zeros(f.N_states, 12);
            G(f.wt, 1:3)   = eye(3)*dt;
            G(f.at, 4:6)   = eye(3)*dt;
            G(f.rbc, 7:9)  = eye(3)*dt;
            G(f.tbc, 10:12)= eye(3)*dt;
            
            noise_nw = 10*eye(3);
            noise_na = 10*eye(3);
            noise_nrc = .1*eye(3);
            noise_ntc = .1*eye(3);
            Q = blkdiag(noise_nw, noise_na, noise_nrc, noise_ntc);  %
            
            P_h = F*P*F' + G*Q*G';
        end
        
        function [x] = states_add(f, dx, x)
           x(f.rsb) = screw_add(dx(f.rsb), x(f.rsb));
           x(f.rbc) = screw_add(dx(f.rbc), x(f.rbc));
           x(f.nonSO3) = dx(f.nonSO3) + x(f.nonSO3); 
        end
    end
end


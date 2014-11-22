classdef ekf_autocalib
    %EKF_AUTOCALIB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X;
        P; % Covariance matrix, same layout as above
        time_stamp;
    end
    
    properties(Constant)
       rsb = 1:3;
       tsb = 4:6;
       wt  = 7:9;
       vt  = 10:12;
%        nt  = 13:15;
       at  = 13:15;
       % Cam related -> should be locked after calibration
       rbc = 16:18; % orientation of camera in body(IMU) frame
       tbc = 19:21; % position    of camera in body(IMU) frame
%        % IMU related
%        wb  = 22:24; % w-bias
%        ab  = 25:27; % a-bias
%        g0  = 28:30; % gravity vector
       nonSO3 = [4:15 19:21];
       N_states = 21;
       Rbc_mean = [1 0 0; 0 -1 0; 0 0 -1];
       max_iter = 10;
       dx_threshold = 1e-16;
    end
    
    methods
        function [f] = ekf_autocalib()
            f.X   = zeros(f.N_states, 1);
            f.P   = 100*eye(f.N_states);

            f.time_stamp = 0.0;
        end
        
        %{
        Process incoming imu data. The measurement model is given as the 
        following:
            y_imu = [    wt    ] + [ wb ] + [ imu_wt_noise ]
                    [R'*(at-g0)]   [ ab ] + [ imu_at_noise ]
            Thus:
            y_pre = [    wt    ] + [ wb ]
                    [R'*(at-g0)]   [ ab ]
        %}
        function [f] = onImuUpdate(f, w, a, imu_noise, time_stamp)
            % Propagate to current time step
            f = propagate(f, time_stamp);
            
            x = f.X; % f.X is the initial guess, stays unchanged for a while
            dx  = zeros(f.N_states, 1);
            dxp = dx;
            H   = zeros(6, f.N_states);
            iter= 0;
            while iter < f.max_iter
                Rsb = screw_exp(x(f.rsb));
                w_pred = x(f.wt) + x(f.wb);
                a_pred = Rsb'*(x(f.at) - x(f.g0)) + x(f.ab);
               
                H(1:3, f.wt) = eye(3);
                H(1:3, f.wb) = eye(3);
                H(4:6, f.at) =  Rsb';
                H(4:6, f.g0) = -Rsb';
                H(4:6, f.rsb)=  Rsb'*so3_alg(f.at-f.g0);
                H(4:6, f.ab) = eye(3);
                
                y = [w - w_pred; a - a_pred] + H*dxp; % iekf innovation
                
                S = H*f.P*H' + imu_noise;
                K = f.P*H'/S;
                dx= K*y;
                
                x = states_add(f.X, dx);                
                if max(abs(dx - dxp)) < f.dx_threshold
                    f.X = x;
                    f.P = f.P - K*H*f.P;
                    converge_flag = TRUE;
                    break;
                end
                dxp = dx;
            end
            
            assert(converge_flag, 'Failed to converge in OnImuUpdate\n');
        end
        
        %{
            y_cam = [  r_sc  ] + [ cam_rx_noise ]
                    [  t_sc  ]   [ cam_cx_noise ]
            R_bc = exp(d_rx_bc)*Rbc_mean;
            where exp(rsc) = exp(rsb)*exp(d_rbc)*Rbc_mean
                      tsc  = exp(rsb)*tbc + tsb
        %}
        
        function [f] = onCamUpdate(f, rsc, tsc, cam_noise, time_stamp)
            % Propagate to current time step
            f = propagate(f, time_stamp);
            
            x = f.X; % f.X is the initial guess, stays unchanged for a while
            dx  = zeros(f.N_states, 1);
            dxp = dx;
            H   = zeros(6, f.N_states);
            iter= 0;
            while iter < f.max_iter
                Rsb = screw_exp(x(f.rsb));
                
                Rsc_pred = Rsb*screw_exp(x(f.rbc))*f.Rbc_mean;
                tsc_pred = Rsb*x(f.tbc) + x(f.tsb);
                
                H(1:3, f.idx_rsb) = eye(3);
                H(1:3, f.idx_rbc) = Rsb;
                H(4:6, f.idx_rsb) = so3_alg(-Rsb*f.tbc);
                H(4:6, f.idx_tbc) = Rsb;
                H(4:6, f.idx_tsb)=  eye(3);
                
                y = [screw_log(screw_exp(rsc)*Rsc_pred');...
                     tsc - tsc_pred] + H*dxp;
                 
                S = H*f.P*H' + cam_noise;
                K = f.P*H'/S;
                dx= K*y;
            
                x = states_add(f.X, dx);
                if max(abs(dx - dxp)) < f.dx_threshold
                    f.X = x;
                    f.P = f.P - K*H*f.P;
                    converge_flag = TRUE;
                    break;
                end
                dxp = dx;
            end
            
            assert(converge_flag, 'Failed to converge in OnCamUpdate\n');
        end
        
        function [f] = propagate(f, new_time)
            
            dt = new_time - f.time_stamp;
            
            % dt could potentially be 0
            if dt < 1e-7
                return;
            end
            
            % TODO: put in a while loop if dt > 1e-1
            
            R_sb = screw_exp(f.X(f.rsb));
            
            f.X(f.rsb) = screw_log(R_sb*screw_exp(f.X(f.wt)*dt));
            f.X(f.tsb) = f.X(f.tsb) + f.X(f.vt) * dt;
            f.X(f.wt)  = f.X(f.wt);%  + f.X(f.nt) * dt;
            f.X(f.vt)  = f.X(f.vt)  + f.X(f.at) * dt;
            
            % Linearized model:
            %{
                 rsb tsb   wt   vt   nt   at  rbc tbc  g0  wb  ab
               F =
                 [I   0   R*dt  0   Rd2t  0    0   0   0   0   0]
                 [0   I    0   I*dt  0   Id2t  0   0   0   0   0]
                 [0   0    I    0   I*dt  0    0   0   0   0   0]
                 [0   0    0    I    0   I*dt  0   0   0   0   0]
                 [0   0    0    0    I    0    0   0   0   0   0]
                 [0   0    0    0    0    I    0   0   0   0   0]
                 [0   0    0    0    0    0    I   0   0   0   0]
                 [0   0    0    0    0    0    0   I   0   0   0]
                 [0   0    0    0    0    0    0   0   I   0   0]
                 [0   0    0    0    0    0    0   0   0   I   0]
                 [0   0    0    0    0    0    0   0   0   0   I]
            Note: Rd2t = R*dt^2/2, Id2t = I*dt^2/2
            %}
            
            F = sparse(eye(f.N_states));
            eye3 = eye(3);
            F(f.rsb, f.wt) = R_sb*dt; %F(f.rsb, f.nt) = R_sb*dt*dt/2;
            F(f.tsb, f.vt) = eye3*dt; %F(f.tsb, f.at) = eye3*dt*dt/2;
            F(f.vt,  f.at) = eye3*dt;
%             F(f.wt,  f.nt) = eye3*dt;
            
            % Uncertainty level proportional to dt:
            wt_noise = 100*dt;
            at_noise = 100*dt;
            g0_noise = 1*dt;
            wb_noise = 1*dt;
            ab_noise = 1*dt;
            rbc_noise= 0*dt;
            tbc_noise= 0*dt;
            Q_prop = sparse(zeros(f.N_states));
%             Q_prop(f.nt, f.nt) = nt_noise*eye3;
            Q_prop(f.wt, f.wt) = wt_noise*eye3;
            Q_prop(f.at, f.at) = at_noise*eye3;
            Q_prop(f.rbc,f.rbc)=rbc_noise*eye3;
            Q_prop(f.tbc,f.tbc)=tbc_noise*eye3;
%             Q_prop(f.g0, f.g0) = g0_noise*eye3;
%             Q_prop(f.wb, f.wb) = wb_noise*eye3;
%             Q_prop(f.ab, f.ab) = ab_noise*eye3;
            
            f.P = F*(f.P + Q_prop)*F';
            f.time_stamp = new_time;
        end

        function [x] = states_add(f, x, dx)
           x(f.rsb) = screw_add(dx(f.rsb), x(f.rsb));
           x(f.rbc) = screw_add(x(f.rbc), dx(f.rbc));
           x(f.nonSO3) = dx(f.nonSO3) + x(f.nonSO3); 
        end
    end     
               
end


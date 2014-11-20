classdef ekf_autocalib
    %EKF_AUTOCALIB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rsb;
        tsb;
        wt;
        vt;
        nt;
        at;
        % IMU related
        g0; % gravity vector
        wb; % w-bias
        ab; % a-bias
        % Cam related -> should be locked after calibration
        rbc; % orientation of camera in body(IMU) frame
        tbc; % position    of camera in body(IMU) frame
        
        P; % Covariance matrix, same layout as above
        time_stamp;
    end
    
    properties(Constant)
       idx_rsb = 1:3;
       idx_tsb = 4:6;
       idx_wt  = 7:9;
       idx_vt  = 10:12;
       idx_nt  = 13:15;
       idx_at  = 16:18;
       idx_g0  = 19:21;
       idx_wb  = 22:24;
       idx_ab  = 25:27;
       idx_rbc = 28:30;
       idx_tbc = 31:33;
    end
    
    methods
        function [f] = ekf_autocalib()
            f.rsb = zeros(3,1);
            f.tsb = zeros(3,1);
            f.wt  = zeros(3,1);
            f.vt  = zeros(3,1);
            f.nt  = zeros(3,1);
            f.at  = zeros(3,1);
            % IMU related
            f.g0  = zeros(3,1); % gravity vector
            f.wb  = zeros(3,1); % w-bias
            f.ab  = zeros(3,1); % a-bias
            % Cam related -> should be locked after calibration
            f.rbc = zeros(3,1); % orientation of camera in body(IMU) frame
            f.tbc = zeros(3,1); % position    of camera in body(IMU) frame
            
            f.P   = 100*eye(33); % Covariance matrix, same layout as above

            f.time_stamp = 0.0;
        end
        
        function [f] = onImuUpdate(f, w, a, imu_noise, time_stamp)
            % Propagate to current time step
            f = propagate(f, time_stamp);
            
            %{
            y_imu = [    wt    ] + [ wb ] + [ imu_wt_noise ]
                    [R'*(at-g0)]   [ ab ] + [ imu_at_noise ]
            Thus:
            y_pre = [    wt    ] + [ wb ]
                    [R'*(at-g0)]   [ ab ]
            %}
            
            R = screw_exp(f.rsb);
            w_pred = f.wt + f.wb;
            a_pred = R'*(f.at - f.g0) + f.ab;
            
            H = zeros(6, 33);
            H(1:3, f.idx_wt) = eye(3);
            H(1:3, f.idx_wb) = eye(3);
            H(4:6, f.idx_at) =  R';
            H(4:6, f.idx_g0) = -R';
            H(4:6, f.idx_rsb)=  R'*so3_alg(f.at-f.g0);
            H(4:6, f.idx_ab) = eye(3);
            
            y = [w - w_pred; a - a_pred];
            S = H*f.P*H' + imu_noise;
            K = f.P*H'/S;
            dx= K*y;
            
            f = add_dx(f, dx);
            f.P = f.P - K*H*f.P;
        end
        
        function [f] = propagate(f, new_time)
            
            dt = new_time - f.time_stamp;
            
            % dt could potentially be 0
            if dt < 1e-7
                return;
            end
            
            % TODO: put in a while loop if dt > 1e-1
            
            R_sb = screw_exp(f.rsb);
            
            f.rsb = screw_log(R_sb*screw_exp(f.wt*dt));
            f.tsb = f.tsb + f.vt * dt;
            f.wt  = f.wt + f.nt * dt;
            f.vt  = f.at + f.at * dt;
            
            % Linearized model:
            %{
                  r   t   wt   vt   nt   at   g0  wb  ab  rbc tbc
               F =
                 [I   0   R*dt  0    0    0    0   0   0   0   0]
                 [0   I    0   I*dt  0    0    0   0   0   0   0]
                 [0   0    I    0   I*dt  0    0   I   0   0   0]
                 [0   0    0    I    0   I*dt  0   0   0   0   0]
                 [0   0    0    0    I    0    0   0   0   0   0]
                 [0   0    0    0    0    I    0   0   I   0   0]
                 [0   0    0    0    0    0    I   0   0   0   0]
                 [0   0    0    0    0    0    0   I   0   0   0]
                 [0   0    0    0    0    0    0   0   I   0   0]
                 [0   0    0    0    0    0    0   0   0   I   0]
                 [0   0    0    0    0    0    0   0   0   0   I]
            %}
            
            F = sparse(eye(33));
            eye3 = eye(3);
            F(f.idx_rsb, f.idx_wt) = R_sb*dt;
            F(f.idx_tsb, f.idx_vt) = eye3*dt;
            F(f.idx_vt,  f.idx_at) = eye3*dt;
            F(f.idx_wt,  f.idx_nt) = eye3*dt;
            F(f.idx_wt,  f.idx_wb) = eye3;
            F(f.idx_at,  f.idx_ab) = eye3;
            
            % Uncertainty level proportional to dt:
            nt_noise = 1*dt;
            at_noise = 1*dt;
            g0_noise = 1*dt;
            wb_noise = 1*dt;
            ab_noise = 1*dt;
            rbc_noise= 1*dt;
            tbc_noise= 1*dt;
            Q_prop = sparse(zeros(33));
            Q_prop(f.idx_nt, f.idx_nt) = nt_noise*eye3;
            Q_prop(f.idx_at, f.idx_at) = at_noise*eye3;
            Q_prop(f.idx_g0, f.idx_g0) = g0_noise*eye3;
            Q_prop(f.idx_wb, f.idx_wb) = wb_noise*eye3;
            Q_prop(f.idx_ab, f.idx_ab) = ab_noise*eye3;
            Q_prop(f.idx_rbc,f.idx_rbc)=rbc_noise*eye3;
            Q_prop(f.idx_tbc,f.idx_tbc)=tbc_noise*eye3;
            
            f.P = F*f.P*F' + Q_prop;
            f.time_stamp = new_time;
        end
        
        function [f] = add_dx(f, dx)
           f.rsb = screw_add(dx(f.idx_rsb), f.rsb);
           f.tsb = dx(f.idx_tsb) + f.tsb;
           f.wt  = dx(f.idx_wt)  + f.wt;
           f.vt  = dx(f.idx_vt)  + f.vt;
           f.nt  = dx(f.idx_nt)  + f.nt;
           f.at  = dx(f.idx_at)  + f.at;
           % IMU related
           f.g0  = dx(f.idx_g0)  + f.g0; 
           f.wb  = dx(f.idx_wb)  + f.wb; 
           f.ab  = dx(f.idx_ab)  + f.ab; 
           % Cam related -> should be locked after calibration
           f.rbc = screw_add(dx(f.idx_rbc), f.rbc); 
           f.tbc = dx(f.idx_tbc) + f.tbc; 
        end
    end
    
end


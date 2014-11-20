classdef autocalibration_naive
    %AUTOCALIBRATION_NAIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X; % [rx cx v w a rbc cbc]
        P;
        time_stamp; % current time in seconds
    end
    
    properties(Constant)
        rx = 1:3;
        tx = 4:6;
        wt = 7:9;
        vt = 10:12;
        nt = 13:15
        at = 16:18;
        rbc = 19:21;  % TODO: change to d_rbc
        tbc = 22:24;
        Rbc_mean = [1  0  0;...
                    0 -1  0;...
                    0  0 -1];
    end
    
    methods
        function f = autocalibration_naive(varargin)
            if nargin == 0
                % Initialize the
                f.X = zeros(24, 0);
                f.P = 100.0*ones(24);
            elseif nargin == 2
                f.X = varargin{1};
                f.P = varargin{2};
            end
            f.time_stamp = 0;
        end
        
        function f = onImuUpdate(f, w, a, R_imu, new_time_stamp)
            % Propagate to current time_stamp:
            f = propagate(f, new_time_stamp);
            
            %{
            y_imu = [   wt   ] + [ imu_wt_noise ]
                    [   at   ]   [ imu_at_noise ]
            but
            y_pre = [   wt   ]
                    [R'*at   ]
            %}
            H = zeros(6, 24);
            Rx = screw_exp(f.X(f.rx));
            a_pred = Rx'*f.X(f.at);
            
            H(1:3, f.wt) = eye(3);
%             H(4:6, f.rx) = so3_alg(a_pred)*Rx'; % see test_diff_so3
            H(4:6, f.at) = Rx';
            
            y = [w - f.X(f.wt);...
                 a -   a_pred];
            S = H*f.P*H' + R_imu;
            K = f.P*H'/S;
            dx= K*y;
            
            f.X(f.rx) = screw_add(dx(f.rx), f.X(f.rx));
            f.X(4:end)= f.X(4:end) + dx(4:end);
            f.P = f.P - K*H*f.P;
        end

        % TODO: Consider iterate this step as in IEKF:
        function f = onCamUpdate(f, r_sc, t_sc, R_cam, time_stamp)
            % Propagate to current time stamp
            f = propagate(f, time_stamp);

            %{
            y_cam = [  rx_sc  ] + [ cam_rx_noise ]
                    [  cx_sc  ]   [ cam_cx_noise ]
            R_bc = exp(d_rx_bc)*Rbc_mean;
            where exp(rx_sc) = exp(rx_sb)*exp(d_rx_bc)*Rbc_mean
                      cx_sc  = exp(rx_sb)*tx_bc + tx_sb
            %}
            R_sb = screw_exp(f.X(f.rx));
            idx_rcb = 1:3;
            idx_tcb = 4:6;
            H = zeros(6, 24);
            H(idx_rcb, f.rx)  = eye(3);
            H(idx_rcb, f.rbc) = R_sb;
            H(idx_tcb, f.rx)  = so3_alg(-R_sb*f.X(f.tbc));
            H(idx_tcb, f.tx)  = eye(3);
            H(idx_tcb, f.tbc) = R_sb;
            
            R_sc_pred = R_sb*screw_exp(f.X(f.rbc))*f.Rbc_mean;
%             r_sc_pred = screw_log(R_sc_pred);
            t_sc_pred = R_sb*f.X(f.tbc) + f.X(f.tx);
            %screw_add(r_sc, -r_sc_pred)
            y = [screw_log(screw_exp(r_sc)*R_sc_pred');...
                 t_sc-t_sc_pred];
            S = H*f.P*H' + R_cam;
            K = f.P*H'/S;
            dx= K*y;
            
            f.X(f.rx) = screw_add(dx(f.rx), f.X(f.rx));
            f.X(4:18) = dx(4:18)  + f.X(4:18); % tx ~ at fields
            f.X(f.rbc)= screw_add(dx(f.rbc),f.X(f.rbc));
            f.X(f.tbc)= dx(f.tbc) + f.X(f.tbc);
            f.P = f.P - K*H*f.P;
        end
        
        function f = propagate(f, new_time_stamp)
            dt = new_time_stamp - f.time_stamp;
            
            % dt could potentially be 0
            if dt < 1e-7
                return;
            end
            
            R = screw_exp(f.X(f.rx));
            f.X(f.rx) = screw_log(R*screw_exp(f.X(f.wt)*dt));
            f.X(f.tx) = f.X(f.tx) + f.X(f.vt)*dt;
            f.X(f.wt) = f.X(f.wt) + f.X(f.nt)*dt;
            f.X(f.vt) = f.X(f.vt) + f.X(f.at)*dt;
            
            % Linearized model:
            %{
                  rx  cx   wt   vt   nt   at   ric tic
               F =
                 [I   0   R*dt  0    0    0    0   0]
                 [0   I    0   I*dt  0    0    0   0]
                 [0   0    I    0   I*dt  0    0   0]
                 [0   0    0    I    0   I*dt  0   0]
                 [0   0    0    0    I    0    0   0]
                 [0   0    0    0    0    I    0   0]
                 [0   0    0    0    0    0    I   0]
                 [0   0    0    0    0    0    0   I]
            %}
            F = eye(24);
            F(f.rx, f.wt) = R*dt;      %F(f.rx, f.nt) = 0.5*R*dt*dt;
            F(f.tx, f.vt) = eye(3)*dt; %F(f.tx, f.at) = 0.5*eye(3)*dt*dt;
            F(f.wt, f.nt) = eye(3)*dt;
            F(f.vt, f.at) = eye(3)*dt;
            
            % Uncertainty level proportional to dt
            angular_noise = 1*dt;
            linear_noise  = 1*dt;
            Q_prop = zeros(24);
            Q_prop(f.nt, f.nt) = angular_noise*eye(3);
            Q_prop(f.at, f.at) =  linear_noise*eye(3);
            
            f.P = F*f.P*F' + Q_prop;
            f.time_stamp = new_time_stamp;
        end
    end
    
end


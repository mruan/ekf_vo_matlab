classdef autocalibration_naive
    %AUTOCALIBRATION_NAIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X; % [rx cx v w a rbi cbi]
        P;
        time_stamp; % current time in seconds
    end
    
    properties(Constant)
        rx = 1:3;
        cx = 4:6;
        wt = 7:9;
        vt = 10:12;
        nt = 13:15
        at = 16:18;
        rbc = 19:21;
        tbc = 22:24;
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
            H(4:6, f.rx) = so3_alg(a_pred)*Rx'; % see test_diff_so3
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
        function f = onCamUpdate(f, rx_sc, cx_sc, R_cam, time_stamp)
            % Propagate to current time stamp
            f = propagate(f, time_stamp);

            %{
            y_cam = [  rx_sc  ] + [ cam_rx_noise ]
                    [  cx_sc  ]   [ cam_cx_noise ]
            where exp(rx_sc) = exp(rx_si)*exp(rx_ic)
                      cx_sc  = exp(rx_si)*cx_ic + cx_si
            %}
            H = zeros(6, 24);
            H(1:3, f.rx) = eye(3);
            H(4:6, f.rx) = so3_alg(-R*cx_ic);
            H(4:6, f.cx) = eye(3);
            
            y = [screw_add(rx_sc, -rx_sc_pred); ...
                           cx_sc - cx_sc_pred];
            S = H*f.P*H' + R_cam;
            K = f.P*H'/S;
            dx= K*y;
            
            f.X(f.rx) = screw_add(dx(f.rx), f.X(f.rx));
            f.X(4:end)= f(4:end) + dx(4:end);
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
            f.X(f.cx) = f.X(f.cx) + f.X(f.vt)*dt;
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
            F(f.rx, f.wt) = R*dt;
            F(f.cx, f.vt) = eye(3)*dt;
            F(f.wt, f.nt) = eye(3)*dt;
            F(f.vt, f.at) = eye(3)*dt;
            
            angular_noise = 1e-2; % 5mm/s^2
            linear_noise  = 1e-2; % 5mm/s^2;
            Q_prop = zeros(24);
            Q_prop(f.nt, f.nt) = angular_noise^2*eye(3);
            Q_prop(f.at, f.at) =  linear_noise^2*eye(3);
            
            f.P = F*f.P*F' + Q_prop;
            f.time_stamp = new_time_stamp;
        end
    end
    
end


function [ ] = w_only_test( )
% Test 1: Can we estimated w from r observations
%   Detailed explanation goes here
dbstop if error

M = gen_data();
t = M.Time;
N = numel(t);

num_state = 9;
X = zeros(num_state, 1);
% X = M.True(1, :)';
P = 100*eye(num_state);

idx_r = 1:3;
idx_w = 4:6;
idx_b = 7:9;

log_X = zeros(num_state, N);
log_P = zeros(num_state, N);

max_iter = 10;
dx_thresh = 1e-4;
t0 = 0;
for i=1:N
    dt = t(i) - t0;
    [X_h, P_h] = propagate(X, P, dt);
    
    if M.Type(i) == 1
        r = M.Data(i,:)';
        r_noise = 1e-8*eye(3);
        [X, P] = Rupdate(r, r_noise, X_h, P_h);
    else
        w = M.Data(i,:)';
        w_noise = 1e-8*eye(3);
        [X, P] = Wupdate(w, w_noise, X_h, P_h);
    end
        
    log_X(:, i) = X - M.True(i,:)';
    log_P(:, i) = sqrt(diag(P));
    t0 = t(i);
    % kalman update
end

subplot(3,1,1);
plot(t, log_X(1,:), 'r', t, log_X(2,:), 'g', t, log_X(3,:), 'b',...
     t, log_P(1,:), 'r-.', t, log_P(2,:), 'g-.', t, log_P(3,:), 'b-.',...
     t,-log_P(1,:), 'r-.', t, log_P(2,:), 'g-.', t,-log_P(3,:), 'b-.');title('r');

subplot(3,1,2);
plot(t, log_X(4,:), 'r', t, log_X(5,:), 'g', t, log_X(6,:), 'b',...
     t, log_P(4,:), 'r-.', t, log_P(5,:), 'g-.', t, log_P(5,:), 'b-.',...
     t,-log_P(4,:), 'r-.', t, log_P(5,:), 'g-.', t,-log_P(5,:), 'b-.');title('w');
 
subplot(3,1,3);
plot(t, log_X(7,:), 'r', t, log_X(8,:), 'g', t, log_X(9,:), 'b',...
     t, log_P(7,:), 'r-.', t, log_P(8,:), 'g-.', t, log_P(9,:), 'b-.',...
     t,-log_P(7,:), 'r-.', t, log_P(8,:), 'g-.', t,-log_P(9,:), 'b-.');title('b');

    function [X, P] = Wupdate(w, noise, X, P)
        x = X; % save the states:
        p = P;
        
        dx = zeros(num_state, 1);
        px = dx;
        H = zeros(3, num_state);
        for iter = 1:max_iter
            w_pred = X(idx_w) + x(idx_b);
            
            H(1:3, idx_w) = eye(3);
            H(1:3, idx_b) = eye(3);
            
            y = w - w_pred;
            
            S = H*P*H' + noise;
            K = P*H'/S;
            dx= K*y;
            
            x = add(dx, x);
            if max(abs(dx - px)) < dx_thresh
                X = x;
%                 P = p - K*S*K';
                P = p - K*H*p;
                break;
            end
            px = dx;
        end% end of loop
    end

    function [X, P] = Rupdate(r, noise, X, P)
        x = X; % save the states:
        p = P;
        
        dx = zeros(num_state, 1);
        px = dx;
        H = zeros(3, num_state);
        for iter = 1:max_iter
            r_pred = x(idx_r);
            
            H(1:3, idx_r) = eye(3);
            
            y = screw_add(r, -r_pred);
            
            S = H*P*H' + noise;
            K = P*H'/S;
            dx= K*y;
            
            x = add(dx, x);
            if max(abs(dx - px)) < dx_thresh
                X = x;
%                 P = p - K*S*K';
                P = p - K*H*p;
                break;
            end
            px = dx;
        end% end of loop
    end

    function [x] = add(dx, x)
       x(idx_r) = screw_add(dx(idx_r), x(idx_r));
       x(4:end) = dx(4:end) + x(4:end);
    end

    function [x, p] = propagate(x, p, dt)
        % propagate
        Rx = screw_exp(x(idx_r));
        w = x(idx_w);
        x(idx_r) = screw_log(Rx*screw_exp(w*dt));
        
        %{
        F = [I 0] + [0 R]
            [0 I]   [0 0]*dt
        %}
        F = eye(num_state);
        F(idx_r, idx_w) = Rx*dt;
        
        Q_w = 10*eye(3);
        Q_b = 0.1*eye(3);
        Q = blkdiag(Q_w, Q_b);
        %{
        G = [0 0]
            [I 0]*dt
            [0 I]
        %}
        G = zeros(num_state, 6);
        G(idx_w, 1:3) = eye(3)*dt;
        G(idx_b, 4:6) = eye(3)*dt;
        
        p = F*p*F' + G*Q*G';
    end
end


function [measurement] = gen_data()
%dt = 0.01;  % 100 Hz
t_end = 10; % 60 seconds
% time_series = 0:dt:t_end;

% Trajecotry control paramters:
measurement.True = [];
measurement.Type = [];
measurement.Time = [];
measurement.Data = [];

wb = [0 0 0]';

IMU_DT = 0.01;
OMN_DT = 0.01;
imu_t = 0.005;
omn_t = 0.00;

while imu_t < t_end || omn_t < t_end
   
    if imu_t <= omn_t
       [a, ad, aD] = time2param(imu_t);
       % IMU measures:
       % y_imu = [    wt    ] + [ wb ] 
       %         [R'*(at-g0)]   [ ab ]
       
       wt = get_wt(a, ad);              
%        at = get_at(a, ad);
       
       Rb = get_R(a);
       rb = screw_log(Rb);
%        Tb = get_T(a);
%        Td = get_Tdot(a, ad);
%        TD = get_Tddot(a, ad, aD);
       
       measurement.Data(end+1, :) = [wt+wb]';
       measurement.True(end+1, :) = [rb; wt; wb]';
       measurement.Type(end+1) = 0;
       measurement.Time(end+1) = imu_t;
       
       imu_t = imu_t + IMU_DT;
    else
        [a, ad, aD] = time2param(omn_t);

        Rb = get_R(a);
%         Tb = get_T(a);
               
        rb = screw_log(Rb);
%         Td = get_Tdot(a, ad);
        wt = get_wt(a, ad);
%         TD = get_Tddot(a, ad, aD);
        
        measurement.Data(end+1, :) = [rb]';
        measurement.True(end+1, :) = [rb; wt; wb]';
        measurement.Type(end+1) = 1;
        measurement.Time(end+1) = omn_t;
        
        omn_t = omn_t + OMN_DT;
    end
    
end
% angle parameter and its time derivatives
    function [a, a_dot, a_ddot] = time2param(t)
        a = t;
        a_dot = 1;
        a_ddot = 0;
        %     omega = pi/10;
        %     C0 = pi/6;
        %     C1 = pi/3;
        %
        %     a = 1.5*omega*t + sin(omega*t+C0) + C1;
        %     a_dot  = omega*(1.5+cos(omega*t+C0));   % always positive
        %     a_ddot = -omega^2*sin(omega*t+C0);
        
        %     a = omega*t;
        %     a_dot = omega*ones(size(t));
        %     a_ddot= zeros(size(t));
    end

function [R] = get_R(a)
    sa = sin(a); ca = cos(a);
    R  = [sa, sa*ca, ca*ca;...
         -ca, sa*sa, ca*sa;...
           0,   -ca,    sa];
end

function [wt] = get_wt(a, a_dot)
    wt = [1; -cos(a); sin(a)]*a_dot;           
end

end
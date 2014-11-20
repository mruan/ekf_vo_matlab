function [full_states, measurement] = gen_env_naive( )
%   
%   Outputs:
%   Trajectory: rx tx, wt, vt, ht, at 
%   Measurement: pose_obs, imu_obs each has an attached timestamp
% Unknowns:
%   1. The relative transformation between IMU and Camera
% Assumptions:
%   1. The pose (rx tx) are directly measureable (implicitly from a camera)
%   2. The IMU readings are compensated,
%      (i.e. correct for drift, exclude gravity and etc)

dt = 0.01;  % 100 Hz
t_end = 60; % 60 seconds
time_series = 0:dt:t_end;

% Trajecotry control paramters:
r = 1;

[full_states] = gen_full_states(time_series);


Rbc = [1  0  0;...
       0 -1  0;...
       0  0 -1]; % rbc -> [pi 0 0];
Tbc = [0.01; 0.03; -0.003];     

[measurement] = gen_measurement(t_end);

save('env_naive.mat', 'full_states', 'measurement');

% angle parameter and its time derivatives
function [a, a_dot, a_ddot] = time2param(t)
    omega = pi/10;
    C0 = pi/6;
    C1 = pi/3;
    a = 1.5*omega*t + sin(omega*t+C0) + C1;
    a_dot  = omega*(1.5+cos(omega*t+C0));   % always positive
    a_ddot = -omega^2*sin(omega*t+C0);

%     a = omega*t;
%     a_dot = omega*ones(size(t));
%     a_ddot= zeros(size(t));
end

% Set up the full states, which capsulate:
% rx: orientation of the body in se(3)
% tx: position of the body
% wt: rotation rate in body frame
% vt: linear velocity in body frame
% at: linear acceleration in body frame
function [full_states] = gen_full_states(t)
num_steps = numel(t);

[a, a_dot, a_ddot] = time2param(t);

full_states.t  = t;
full_states.a  = a;
full_states.ad = a_dot;
full_states.aD = a_ddot;
full_states.tx = get_T (r, a);
full_states.vt = get_vt(r, a, a_dot); 
full_states.at = get_at(r, a, a_dot, a_ddot);
full_states.wt = get_wt(a, a_dot);

full_states.rx = zeros(3, num_steps);
for i=1:num_steps
    R = get_R(a(i));
    full_states.rx(:, i) = screw_log(R);
end % end of for

end % end of FUNCTION:gen_full_states()

% Measurement has 3 fields: TimeStamp, Type and Data
% TimeStamp: Time the observation becomes available
% Type(int): 0 for imu, 1 for camera
% Data(1x6): [w a]' for imu, [rx tx]' for camera
% Cheat(15x1): [rx(1-3) tx(4-6) wt(7-9) vt(10-12) at(13-15)]
function [measurement] = gen_measurement(t_end)
    IMU_DT = 0.01; % 100 Hz
    CAM_DT = 1/30; %  30 Hz
    
    measurement.Type = [];
    measurement.TimeStamp = [];
    measurement.Data = [];
    measurement.Cheat = [];
    
    imu_t = IMU_DT; % slight offset
    cam_t = 0.0;
    
    while imu_t < t_end && cam_t < t_end
        
        % it possible that both sensors update at the same time
        if imu_t <= cam_t            
            [a, a_dot, a_ddot] = time2param(imu_t);
            
            % Used in data
            wt = get_wt(a, a_dot);
            at = get_at(r, a, a_dot, a_ddot);
            
            % Used in cheats
            Rb = get_R(a); 
            rb = screw_log(Rb);
            Tb = get_T(r, a);
            Td = get_Tdot(r, a, a_dot);
            Tdd= Rb*at;
            
            measurement.Type(end+1)      = 0;
            measurement.TimeStamp(end+1) = imu_t;
            measurement.Data(end+1,:)    = [wt; at]';
            measurement.Cheat(end+1,:)   = [rb; Tb; wt; Td; Tdd]';

            imu_t = imu_t + IMU_DT;
        else
            
            [a, a_dot, a_ddot] = time2param(cam_t);
            
            % Used in data
            Rb = get_R(a);
            Tb = get_T(r, a);
            
            Rc = Rb*Rbc;        % Camera orientation
            rc = screw_log(Rc);
            Tc = Rb*Tbc + Tb;   % Camera position

            % Used in cheat
            rb = screw_log(Rb);
            Td = get_Tdot(r, a, a_dot);
            wt = get_wt(a, a_dot);
            Tdd= Rb*get_at(r, a, a_dot, a_ddot);
            
            measurement.Type(end+1)      = 1;
            measurement.TimeStamp(end+1) = cam_t;
            measurement.Data(end+1,:)    = [rc; Tc]';            
            measurement.Cheat(end+1,:)   = [rb; Tb; wt; Td; Tdd]';
            
            cam_t = cam_t + CAM_DT;
        end
    end
    
end % end of FUNCTION:gen_measurement()

function [T] = get_T(r, a)
    T = [cos(a);...
            0*a;...
         sin(a)]*r;
end

function [Td] = get_Tdot(r, a, a_dot)
    Td = [-sin(a); ...
              0*a; ...
           cos(a)]*r.*a_dot;
end

function [R] = get_R(a)
    sa = sin(a); ca = cos(a);
    R  = [ca/sqrt(2), -sa, ca/sqrt(2);...
           1/sqrt(2),   0, -1/sqrt(2);...
          sa/sqrt(2),  ca, sa/sqrt(2)];
end

function [vt] = get_vt(r, a, a_dot)
    vt = [    0*a; ...
          r*a_dot; ...
              0*a];
end

function [wt] = get_wt(a, a_dot)
    wt = [-1*a_dot;...
           0*a_dot;...
           1*a_dot]/sqrt(2);            
end

function [at] = get_at(r, a, a_dot, a_ddot)
    at = [-a_dot.*a_dot/sqrt(2);...
                         a_ddot;...
          -a_dot.*a_dot/sqrt(2)]*r;
end

end % end of CODE

%{
Let's assume the sensor moves in a circular motion
   centered at [0 0 0], with radius r, in the x-z plane 
Orientation-wise, it looks at a fixed point [0 R 0]
   the motion parameter, a, is the angle of the cicular motion, and is
   a function of t (so that acceleration is non-zero)

a = sym('a', [1 1]);
a_dot = sym('a_dot', [1 1]);
a_ddot= sym('a_ddot', [1 1]);
r = sym('r', [1 1]);

Body Pose states:
R = [
[ cos(a)/sqrt(2), -sin(a),  cos(a)/sqrt(2)]
[      1/sqrt(2),       0,      -1/sqrt(2)]
[ sin(a)/sqrt(2),  cos(a),  sin(a)/sqrt(2)]]

T = [
[ cos(a)]
[     0 ]
[ sin(a)]] * r

First order derivative:
R_dot = [                                              
[ -sin(a)/sqrt(2), -cos(a), -sin(a)/sqrt(2)]           
[               0,       0,               0]   
[  cos(a)/sqrt(2), -sin(a),  cos(a)/sqrt(2)]] * a_dot     

% T_dot =def= (R*v_sb^b)
T_dot = [
[-sin(a)]
[     0 ]
[ cos(a)]] * r * a_dot

% Velocity states:
w_sb_b_x = simplify(transpose(R)*R_dot) 
% [w_sb^b]_x = R'*R_dot
% [w_sb^b]_x = [           
% [0, -1,  0]                           
% [1,  0,  1]
% [0, -1,  0]]/sqrt(2)*a_dot

w_sb_b = so3_gla(w_sb_b_x) 
%[-1]
%[ 0]/sqrt(2)*a_dot;
%[ 1]

v_sb_b = transpose(R) * T_dot
%[ 0 ]
%[ 1 ]*r*a_dot
%[ 0 ]

Acceleration states: 
Definition: v_sb_b_dot = diff(v_sb_b)
v_sb_b_dot = [0; 1; 0]*r*a_ddot
Also known:
v_sb_b_dot = -[w_sb^b]_x * v_sb^b + a_sb^b
a_sb_b = v_sb_b_dot + w_sb_b_x * v_sb_b=
[ 0]             [-1]
[ 1]*r*a_ddot +  [ 0]/sqrt(2)*r*a_dot^2
[ 0]             [-1]
(in terminal:)
a_sb_b = [
[-a_dot.*a_dot/sqrt(2);...
                a_ddot;...
 -a_dot.*a_dot/sqrt(2)]]*r

Measurements:
y_imu = [w_sb^b; a_sb^b] + N(0, Q_imu)
y_cam = [R T]_sb*[R T]_bc
%}




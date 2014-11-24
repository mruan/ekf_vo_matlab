function [measurement] = gen_sim_camcalib( )
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

dbstop if error

%dt = 0.01;  % 100 Hz
t_end = 10; % 60 seconds
% time_series = 0:dt:t_end;

% Trajecotry control paramters:
r = 1;

% Camera pose in body frame
Rbc = [1  0  0;...
       0 -1  0;...
       0  0 -1]; % rbc -> [pi 0 0];
rbc = screw_log(Rbc);
tbc = [0.01; 0.03; -0.003];     

% Measurement has 4 fields: Time, Type True and Data
% Time(double): Time the observation becomes available
% Type(int)   : 0 for imu, 1 for camera
% Data(1x6)   : [w a]' for imu, [rx tx]' for camera
% True(15x1)  : [rx(1-3) tx(4-6) wt(7-9) Td(10-12) Tdd(13-15)]

IMU_DT = 0.01;
CAM_DT = 0.01;

measurement.True = [];
measurement.Type = [];
measurement.Time = [];
measurement.Data = [];

imu_t = 0.0;
cam_t = 0.005;

while imu_t < t_end || cam_t < t_end
   
    if imu_t <= cam_t
       [a, ad, aD] = time2param(imu_t);
       % IMU measures:
       % y_imu = [    wt    ] + [ wb ] 
       %         [R'*(at-g0)]   [ ab ]
       
       wt = get_wt(a, ad);              
       at = get_at(a, ad);
       
       Rb = get_R(a);
       rb = screw_log(Rb);
       Tb = get_T(a);
       Td = get_Tdot(a, ad);
       TD = get_Tddot(a, ad, aD);
       
       measurement.Data(end+1, :) = [wt; at]'; %
       measurement.True(end+1, :) = [rb; Tb; wt; Td; TD; rbc; tbc];
       measurement.Type(end+1) = 0;
       measurement.Time(end+1) = imu_t;
       
       imu_t = imu_t + IMU_DT;
    else
        [a, ad, aD] = time2param(cam_t);

        Rb = get_R(a);
        Tb = get_T(a);
        
        Rsc = Rb*Rbc;
        rsc = screw_log(Rsc);
        tsc = Rb*tbc + Tb;
        
        rb = screw_log(Rb);
        Td = get_Tdot(a, ad);
        wt = get_wt(a, ad);
        TD = get_Tddot(a, ad, aD);
        
        measurement.Data(end+1, :) = [rsc; tsc]';%[rb; Tb]';
        measurement.True(end+1, :) = [rb; Tb; wt; Td; TD; rbc; tbc];
        measurement.Type(end+1) = 1;
        measurement.Time(end+1) = cam_t;
        
        cam_t = cam_t + CAM_DT;
    end
    
end

save('sim_env_cam.mat', 'measurement');

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

% % Set up the full states, which capsulate:
% % rx: orientation of the body in se(3)
% % tx: position of the body
% % wt: rotation rate in body frame
% % vt: linear velocity in body frame
% % at: linear acceleration in body frame
% function [full_states] = gen_full_states(t)
% num_steps = numel(t);
% 
% [a, a_dot, a_ddot] = time2param(t);
% 
% full_states.t  = t;
% full_states.a  = a;
% full_states.ad = a_dot;
% full_states.aD = a_ddot;
% full_states.tx = get_T (r, a);
% full_states.vt = get_vt(r, a, a_dot); 
% full_states.at = get_at(r, a, a_dot, a_ddot);
% full_states.wt = get_wt(a, a_dot);
% 
% full_states.rx = zeros(3, num_steps);
% for i=1:num_steps
%     R = get_R(a(i));
%     full_states.rx(:, i) = screw_log(R);
% end % end of for
% 
% end % end of FUNCTION:gen_full_states()

function [R] = get_R(a)
    sa = sin(a); ca = cos(a);
    R  = [sa, sa*ca, ca*ca;...
         -ca, sa*sa, ca*sa;...
           0,   -ca,    sa];
end

function [T] = get_T(a)
    T = [0 0 0]'*a;
end

% Note: THIS IS NOT SAME AS v_sb^b
function [Td] = get_Tdot(a, a_dot)
    Td = [0 0 0]'*a;
end

% NOTE: THIS IS NOT SAME AS a_sb^b
function [Tdd] = get_Tddot(a, a_dot, a_ddot)
    Tdd= [0 0 0]'*a;
end

% function [vt] = get_vt(r, a, a_dot)
%     vt = [    0*a; ...
%           r*a_dot; ...
%               0*a];
% end
% 
function [wt] = get_wt(a, a_dot)
    wt = [1; -cos(a); sin(a)]*a_dot;           
end

function [at] = get_at(a, a_dot, a_ddot)
    at = [0 0 0]'*a;
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

Body Pose states:
R = [
[ sin(a), sin(a).*cos(a),  cos(a).*cos(a)]
[-cos(a), sin(a).*sin(a),  cos(a).*sin(a)]
[     0 ,        -cos(a),          sin(a)]]
simplify(det(R)) (==1)
simplify(transpose(R)*R) ==eye(3)

First order derivative:
R_dot = [                                              
[ cos(a), cos(2*a), -sin(2*a)]           
[ sin(a), sin(2*a),  cos(2*a)]   
[     0 ,   sin(a),    cos(a)]]*a_dot

=> w_sb_b_x = simplify(transpose(R)*R_dot)
=> w_sb_b   = so3_gla(w_sb_b_x) 
=> w_sb_b   = [1; -cos(a); sin(a)]*a_dot

T = [0 0 0];
Td= [0 0 0];
TD= [0 0 0];

From Td = def= (R*v_sb^b)
v_sb_b = transpose(R)*Td = [0 0 0];
v_sb_b_dot = diff(v_sb_b) => [0 0 0];
a_sb_b     = v_sb_b_dot + w_sb_b_x * v_sb_b = [0 0 0]

Measurements:
y_imu = [w_sb^b; a_sb^b] + N(0, Q_imu)
y_cam = [R T]_sb*[R T]_bc
%}
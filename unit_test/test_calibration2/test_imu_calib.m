% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

%clear;
function [f] = test_imu_calib()
load('sim_env2.mat'); % measurement
addpath('../../lie_group');

t = measurement.Time;
N = numel(t);

rx = zeros(3, N);
tx = zeros(3, N);
vt = zeros(3, N);
at = zeros(3, N);
rbc= zeros(3, N);
tbc= zeros(3, N);

f = ekf_autocalib();
f.X(f.rsb) = measurement.True(1, 1:3)'; f.P(f.rsb, f.rsb) = 1e-2*eye(3);
f.X(f.tsb) = measurement.True(1, 4:6)'; f.P(f.tsb, f.tsb) = 1e-2*eye(3);
% f.X(f.wt)  = measurement.True(1, 7:9)';   f.P(f.wt, f.wt) = 1e-8*eye(3);
% f.X(f.vt)  = measurement.True(1, 10:12)'; f.P(f.vt, f.vt) = 1e-8*eye(3);
% f.X(f.at)  = measurement.True(1, 13:15)'; f.P(f.at, f.at) = 1e-8*eye(3);
f.X(f.rbc) = [0.05 0.05 0.05]'; f.P(f.rbc, f.rbc) = 9e-4*eye(3);
f.X(f.tbc) = [0.00; 0.00; -0.000]; f.P(f.tbc, f.tbc) = 9e-4*eye(3);

% f.X(f.g0 ) = [0 0 9.8]';
%f.wt  = measurement.True(1, 7:9)';
%f.vt  = measurement.True(1, 10:12)';

for i=1:N
    
    comp = @(x, i) [x measurement.True(i,:)' x-measurement.True(i,:)' [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7]'];
    if measurement.Type(i) == 0
        w = measurement.Data(i, 1:3)';
        a = measurement.Data(i, 4:6)';
        
        imu_noise = 1e-6*eye(6);

        current_time = t(i);
        f = propagate(f, current_time);
        
            x = f.X; % f.X is the initial guess, stays unchanged for a while
            dx  = zeros(f.N_states, 1);
            comp_x = comp(x, i); %IMU
            dxp = dx;
            H   = zeros(6, f.N_states);
            converge_flag = 0;
            for iter = 1:f.max_iter
                Rsb = screw_exp(x(f.rsb));
                w_pred = x(f.wt);
                a_pred = Rsb'*x(f.at);
               
                H(1:3, f.wt) = eye(3);
                H(4:6, f.at) =  Rsb';
                H(4:6, f.rsb)=  Rsb'*so3_alg(x(f.at));
                
                y = [w - w_pred; a - a_pred] + H*dxp; % iekf innovation
                
                S = H*f.P*H' + imu_noise;
                K = f.P*H'/S;
                dx= K*y;
                
                x = f.states_add(f.X, dx);                
                diff = dx - dxp;
                if max(abs(diff)) < f.dx_threshold
                    f.X = x;
                    f.P = f.P - K*H*f.P;
                    converge_flag = 1;
%                     comp_x = [extract_sub(x) measurement.True(i,:)' [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7]'];
             comp_x = comp(x, i); %IMU
                    break;
                end
                dxp = dx;
            end
            
            
%            assert(converge_flag>0, 'Failed to converge in OnImuUpdate\n');
    else % Camera update
       
        rsc = measurement.Data(i, 1:3)';
        tsc = measurement.Data(i, 4:6)';
        
        cam_noise = 1e-6*eye(6);

            current_time = t(i);
            f = propagate(f, current_time);

            x = f.X; % f.X is the initial guess, stays unchanged for a while
            comp_x = comp(x, i); % CAM
            dx  = zeros(f.N_states, 1);
            dxp = dx;
            H   = zeros(6, f.N_states);
            converge_flag = 0;
            for iter = 1:f.max_iter
                Rsb = screw_exp(x(f.rsb));
                
                Rsc_pred = Rsb*screw_exp(x(f.rbc))*f.Rbc_mean;
                tsc_pred = Rsb*x(f.tbc) + x(f.tsb);
                
                H(1:3, f.rsb) = eye(3);
                H(1:3, f.rbc) = Rsb;
                H(4:6, f.rsb) = -so3_alg(Rsb*x(f.tbc));
                H(4:6, f.tbc) = Rsb;
                H(4:6, f.tsb)=  eye(3);
                
                y = [screw_log(screw_exp(rsc)*Rsc_pred');...
                     tsc - tsc_pred] + H*dxp;
                 
                S = H*f.P*H' + cam_noise;
                K = f.P*H'/S;
                dx= K*y;
            
                x = f.states_add(f.X, dx);
                diff = dx - dxp;
                if max(abs(diff)) < f.dx_threshold
                    f.X = x;
                    f.P = f.P - K*H*f.P;
                    converge_flag = 1;
              comp_x = comp(x, i); % CAM
                    break;
                end
                dxp = dx;
            end
            
%            assert(converge_flag>0, 'Failed to converge in OnCamUpdate\n');
    end
    
    rx(:,i) = f.X(f.rsb);
    tx(:,i) = f.X(f.tsb);
    vt(:,i) = f.X(f.vt);
    wt(:,i) = f.X(f.wt);
    at(:,i) = f.X(f.at);
    rbc(:, i) = f.X(f.rbc);
    tbc(:, i) = f.X(f.tbc);
end

subplot(3,2,1);
plot(t, rx(1,:)'-measurement.True(:, 1), 'r', ...
     t, rx(2,:)'-measurement.True(:, 2), 'g', ...
     t, rx(3,:)'-measurement.True(:, 3), 'b'); title('rx');
 
subplot(3,2,2);
plot(t, tx(1,:)'-measurement.True(:, 4), 'r', ...
     t, tx(2,:)'-measurement.True(:, 5), 'g', ...
     t, tx(3,:)'-measurement.True(:, 6), 'b'); title('tx');
 
subplot(3,2,3);
plot(t, vt(1,:)'-measurement.True(:, 10), 'r', ...
     t, vt(2,:)'-measurement.True(:, 11), 'g', ...
     t, vt(3,:)'-measurement.True(:, 12), 'b'); title('vt');
subplot(3,2,4);
plot(t, at(1,:)'-measurement.True(:, 13), 'r', ...
     t, at(2,:)'-measurement.True(:, 14), 'g', ...
     t, at(3,:)'-measurement.True(:, 15), 'b'); title('at'); 
 
subplot(3,2,5);
% plot(t, vt(1,:)'-measurement.True(:, 10), 'r', ...
%      t, vt(2,:)'-measurement.True(:, 11), 'g', ...
%      t, vt(3,:)'-measurement.True(:, 12), 'b'); title('vt');
plot(t, rbc(1, :), 'r', t, rbc(2,:), 'g', t, rbc(3,:)); title('rbc');

subplot(3,2,6);
% plot(t, wt(1,:)'-measurement.True(:, 7), 'r', ...
%      t, wt(2,:)'-measurement.True(:, 8), 'g', ...
%      t, wt(3,:)'-measurement.True(:, 9), 'b'); title('wt');
plot(t, tbc(1, :), 'r', t, tbc(2,:), 'g', t, tbc(3,:)); title('tbc');
end

function [x] = extract_sub(x)
   x = x([1:12 16:18 19:24]);
end
% subplot(3,1,1);
% plot(t, measurement.True(:, 1), 'r', ...
%      t, measurement.True(:, 2), 'g', ...
%      t, measurement.True(:, 3), 'b'); hold on;
% plot(t, rx(1,:), 'r-.', t, rx(2,:), 'g-.', t, rx(3,:), 'b-.'); hold off;
% 
% subplot(3,1,2);
% plot(t, measurement.True(:, 4), 'r', ...
%      t, measurement.True(:, 5), 'g', ...
%      t, measurement.True(:, 6), 'b'); hold on;
% plot(t, tx(1,:), 'r-.', t, tx(2,:), 'g-.', t, tx(3,:), 'b-.'); hold off;
% 
% subplot(3,1,3);
% plot(t, measurement.True(:, 10), 'r', ...
%      t, measurement.True(:, 11), 'g', ...
%      t, measurement.True(:, 12), 'b'); hold on;
% plot(t, vt(1,:), 'r-.', t, vt(2,:), 'g-.', t, vt(3,:), 'b-.'); hold off;

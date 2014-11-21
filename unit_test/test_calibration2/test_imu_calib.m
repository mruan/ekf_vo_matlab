% Load the generated sim_env, then verify if True is indeed ground truth

% Given ground truth at, and using constant acceleration model,
% position drift with time.

%clear;
function [f] = test_imu_calib()
load('sim_env.mat'); % measurement

t = measurement.Time;
N = numel(t);

rx = zeros(3, N);
tx = zeros(3, N);
vt = zeros(3, N);
at = zeros(3, N);
rbc= zeros(3, N);
tbc= zeros(3, N);

f = ekf_autocalib();
f.X(f.rsb) = measurement.True(1, 1:3)';
f.X(f.tsb) = measurement.True(1, 4:6)';
f.X(f.g0 ) = [0 0 9.8]';
%f.wt  = measurement.True(1, 7:9)';
%f.vt  = measurement.True(1, 10:12)';

for i=1:N
    
    if measurement.Type(i) == 0
        w = measurement.Data(i, 1:3)';
        a = measurement.Data(i, 4:6)';
        
        imu_noise = 1e-6*eye(6);
%         
%         f = propagate(f, t(i));

% Test what if the filter has omni-knowledge of the states:
%         f.X(f.rsb) = measurement.True(i, 1:3)';   f.P(f.rsb, f.rsb) = 1e-6*eye(3);
%         f.X(f.tsb) = measurement.True(i, 4:6)';   f.P(f.tsb, f.tsb) = 1e-6*eye(3);
        f.X(f.wt)  = measurement.True(i, 7:9)';   f.P(f.wt, f.wt) = 1e-8*eye(3);
        f.X(f.vt)  = measurement.True(i, 10:12)'; f.P(f.vt, f.vt) = 1e-8*eye(3);
        f.X(f.at)  = measurement.True(i, 13:15)'; f.P(f.at, f.at) = 1e-8*eye(3);
    
            x = f.X; % f.X is the initial guess, stays unchanged for a while
            dx  = zeros(f.N_states, 1);
            dxp = dx;
            H   = zeros(6, f.N_states);
            converge_flag = 0;
            for iter = 1:f.max_iter
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
                
                x = f.states_add(f.X, dx);                
                diff = dx - dxp;
                if max(abs(diff)) < f.dx_threshold
                    f.X = x;
                    f.P = f.P - K*H*f.P;
                    converge_flag = 1;
                    break;
                end
                dxp = dx;
            end
            
            estm_x = extract_sub(x);
            true_x = measurement.True(i,:)';
            assert(converge_flag>0, 'Failed to converge in OnImuUpdate\n');
    else % Camera update
       
        rsc = measurement.Data(i, 1:3)';
        tsc = measurement.Data(i, 4:6)';
        
        cam_noise = 1e-6*eye(6);


%         f.rsb = measurement.True(i, 1:3)';   f.P(f.idx_rsb, f.idx_rsb) = 1e-6*eye(3);
%         f.tsb = measurement.True(i, 4:6)';   f.P(f.idx_tsb, f.idx_tsb) = 1e-6*eye(3);
%         f.X(f.wt)  = measurement.True(i, 7:9)';   f.P(f.wt, f.wt) = 1e-8*eye(3);
%         f.X(f.vt)  = measurement.True(i, 10:12)'; f.P(f.vt, f.vt) = 1e-8*eye(3);
%         f.X(f.at)  = measurement.True(i, 13:15)'; f.P(f.at, f.at) = 1e-8*eye(3);
            f = propagate(f, t(i));

            x = f.X; % f.X is the initial guess, stays unchanged for a while
            estm_x = extract_sub(x);
            true_x = measurement.True(i,:)';
            
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
                H(4:6, f.rsb) = so3_alg(-Rsb*x(f.tbc));
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
                    break;
                end
                dxp = dx;
            end
            
            assert(converge_flag>0, 'Failed to converge in OnCamUpdate\n');
        
        r_bc = f.X(f.rbc);
        t_bc = f.X(f.tbc);
%         f.onCamUpdate;(rsc, tsc, cam_noise, t(i));
    end
    
    rx(:,i) = f.rsb;
    tx(:,i) = f.tsb;
    vt(:,i) = f.vt;
    wt(:,i) = f.wt;
    at(:,i) = f.at;
    rbc(:, i) = f.rbc;
    tbc(:, i) = f.tbc;
end

subplot(2,2,1);
plot(t, rx(1,:)'-measurement.True(:, 1), 'r', ...
     t, rx(2,:)'-measurement.True(:, 2), 'g', ...
     t, rx(3,:)'-measurement.True(:, 3), 'b'); title('rx');
 
subplot(2,2,2);
plot(t, tx(1,:)'-measurement.True(:, 4), 'r', ...
     t, tx(2,:)'-measurement.True(:, 5), 'g', ...
     t, tx(3,:)'-measurement.True(:, 6), 'b'); title('tx');
 
subplot(2,2,3);
% plot(t, vt(1,:)'-measurement.True(:, 10), 'r', ...
%      t, vt(2,:)'-measurement.True(:, 11), 'g', ...
%      t, vt(3,:)'-measurement.True(:, 12), 'b'); title('vt');
plot(t, rbc(1, :), 'r', t, rbc(2,:), 'g', t, rbc(3,:)); title('rbc');

subplot(2,2,4);
% plot(t, wt(1,:)'-measurement.True(:, 7), 'r', ...
%      t, wt(2,:)'-measurement.True(:, 8), 'g', ...
%      t, wt(3,:)'-measurement.True(:, 9), 'b'); title('wt');
plot(t, tbc(1, :), 'r', t, tbc(2,:), 'g', t, tbc(3,:)); title('tbc');
end

function [x] = extract_sub(x)
   x = x([1:12 16:18]);
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

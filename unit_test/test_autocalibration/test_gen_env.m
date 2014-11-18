% Test gen_env:
% Numerically integrate the states to recover path and measurements
% If they don't match reasonably well, that's gonna suck...

load('env_naive.mat');

% the time scale
t  = full_states.t;
dt = t(2) - t(1);
num_steps = numel(t);

% Numerically integrate states:
rx = zeros(3, num_steps);
tx = zeros(3, num_steps);
vt = zeros(3, num_steps);
wt = zeros(3, num_steps);
at = zeros(3, num_steps);

rx(:,1) = full_states.rx(:,1);
tx(:,1) = full_states.tx(:,1);
R_im1 = screw_exp(rx(:,1));
vt(:,1) = full_states.vt(:,1);
wt(:,1) = full_states.wt(:,1);
at(:,1) = full_states.at(:,1);

for i= 2:num_steps   
    R_im  = R_im1*screw_exp(wt(:,i-1)*dt);
    
    rx(:,i) = screw_log(R_im);
    tx(:,i) = tx(:,i-1) + R_im1*vt(:,i-1)*dt;
    
    dv = -so3_alg(wt(:,i-1))*vt(:,i-1)+at(:,i-1);
    vt(:,i) = vt(:,i) + dv*dt;
    
    wt(:,i) = full_states.wt(:,i);
    at(:,i) = full_states.at(:,i);
    
    R_im1 = R_im;
    
%     R_im  = R_im1*screw_exp(wt(:,i-1)*dt);
%     
%     rx(:,i) = screw_log(R_im);
%     tx(:,i) = tx(:,i-1) + vt(:,i-1)*dt;
% 
%     vt(:,i) = vt(:,i-1) + at(:,i-1)*dt;
%     
%     wt(:,i) = full_states.wt(:,i);
%     
%     at(:,i) = R_im*full_states.at(:,i);
%     R_im1 = R_im;    
end
% Plot the states:
plot(1); clf;
subplot(2,1,1);
plot(t, rx(1,:), 'r', t, rx(2,:), 'g', t, rx(3,:), 'b',...
     t, full_states.rx(1,:), 'r-.', ...
     t, full_states.rx(2,:), 'g-.', ...
     t, full_states.rx(3,:), 'b-.');

subplot(2,1,2);
plot(t, tx(1,:), 'r', t, tx(2,:), 'g', t, tx(3,:), 'b',...
     t, full_states.tx(1,:), 'r-.', ...
     t, full_states.tx(2,:), 'g-.', ...
     t, full_states.tx(3,:), 'b-.');
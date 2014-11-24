
%% Seed generator
rng(1);

%% Test exp(r_t') = exp(r_t)exp(w*dt)
r_t = rand(3,1);
w   = rand(3,1);
dt  = 1e-5;

R_tp = expm(so3_alg(r_t))*expm(so3_alg(w*dt));
r_tp = so3_gla(logm(R_tp));
delta = 1e-5;
for i=1:3
    di = zeros(3, 1);
    di(i) = delta;
    
    R_tp_del = expm(so3_alg(di))*expm(so3_alg(r_t))*expm(so3_alg(w*dt));
    
    do = so3_gla(logm(R_tp_del/R_tp));
    Hr(:,i) = do/delta;
end

Hr

for i=1:3
    di = zeros(3, 1);
    di(i) = delta;
    
    R_tp_del = expm(so3_alg(r_t))*expm(so3_alg(di*dt))*expm(so3_alg(w*dt));
    
    do = so3_gla(logm(R_tp_del/R_tp));
    Hw(:,i) = do/delta;
end
Hw
expm(so3_alg(r_t))*dt

%% Test the following derivative:
%   d/dx exp(x)'*a
%  =d/dx exp(-x)*a
%  =d/dx exp(-(z<+>x))
%  =d/dz exp(-x)*exp(-z)*a
% ...
%  = [exp(x)'*a]_x*exp(x)'
% rx = rand(3,1);
% a = rand(3,1);
% Rx = screw_exp(rx);
% 
% delta = 1e-6;
% Hn = zeros(3,3);
% 
% for i=1:3
%    d_rx = [0 0 0]';
%    d_rx(i) = delta;
%    
%    Rx_p = screw_exp(d_rx)*Rx;
%    
%    ai = Rx'  *a;
%    ao = Rx_p'*a;
%    
%    Hn(:,i) = (ao-ai)/delta;
% end
% 
% Hn
% 
% so3_alg(Rx'*a)*Rx'
% Rx'*so3_alg(a)


% %% Test continuous time derivative:
% rx_t0 = rand(3,1);
% w  = rand(3,1);
% dt = 1e-2;
% Rx = screw_exp(rx);
% Drx = Rx*screw_exp(w*dt);
% 
% rx_t1 = screw_add(rx_t0, w*dt);
% 
% d_rx = screw_add(rx_t1, -rx_t0)/dt
% 
% d_rx_a = screw_exp(rx_t0)*w

% delta = 1e-6;
% Hn= zeros(3,3);
% for i=1:3
%     d_rx = [0 0 0]';
%     d_rx(i) = delta;
%     
%     Rx_p = screw_exp(d_rx)*Rx;
%     Drxp = Rx_p*so3_alg(w);
%     
%     Hn(:,i) = screw_add(Drxp, -Drx)/delta;
% end
% 
% Hn
% so3_alg(w)



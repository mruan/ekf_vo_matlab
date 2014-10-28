
rng(0);

%% Generate some random twist from a Gaussian:
theta_1 = 10*rand(3,1);

theta_2 = 10*rand(3,1);

R1 = screw_exp(theta_1);

R2 = screw_exp(theta_2);

R3 = R2*R1;

theta_3 = screw_log(R3); 

H1 = zeros(3,3);
delta = 1e-6;
for i=1:3
    theta = zeros(3,1);
    theta(i) = delta;
    R1p = screw_exp(theta)*R1;
    
    R3p = R2*R1p;
    
    delta_theta = screw_log(R3p*R3');
    
    H1(:, i) = delta_theta/delta;
end

H1
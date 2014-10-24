

K = [500 0 320; 0 500 240; 0 0 1];
%% Point at inifinity:
p_inf = [0 30 200]';

target = p_inf;
%% Horizontal parallex:
N = 11;
x = linspace(-5, 5, N);
y = zeros(N,1);
z = -5*ones(N,1);
M = zeros(3, 4, N);
for i=1:N
   M(:,1:3,i) = eye(3);
   M(:,  4,i) = -M(:,1:3,i)*[x(i) y(i) z(i)]';
end

%% We want random but repeatable noise
rng(1);

filter = ekf_P3;
for i=1:N
   %% Project the world point to image 
   P = M(:,1:3, i)*target+M(:,4,i);
   p = K*P;
    
   z = [p(1)/p(3) p(2)/p(3)]';
   
   %% (+additive isometric Gaussian noise)
   z = z + 0.33*randn(2,1);
   
   %% Initialize a new point or run filter
   if i==1
       % Initialize the Kalman Filter with the first image
       filter = filter.init(  z, K, M(:,1:3,i), M(:,4,i));
   else
       % Run the filter
       filter = filter.update(z, K, M(:,1:3,i), M(:,4,i));
   end
   fprintf('Iteration %d, info:\n', i);
   filter.x
   filter.P
end
    
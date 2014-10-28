
%% Part 1: Test vector and skew-symmetric matrix conversion
w = rand(3,1);
% Test vector to skew-symmetric matrix
m = so3_alg(w)
m + m' % this must add up to zeros

% Test skew-symmetric matrix to vector
wp = so3_gla(m)
w

%% Part 2: Test screw exp and log:
fprintf('Screw: Rodriguez test');
w = rand(3,1);
w = 1.5*w/norm(w);
expm_Result = expm(so3_alg(w))
exp_Result  = screw_exp(w)
wp= screw_log(exp_Result)
w

fprintf('Taylor test');
w = 0.1*w/norm(w);
expm_Result = expm(so3_alg(w))
exp_Result  = screw_exp(w)
wp= screw_log(exp_Result)
w

%% Part 3: Test twist exp and log:
fprintf('Twist: Rodriguez test');
w = rand(6,1);
expm_Result = expm(se3_alg(w))
exp_Result  = twist_exp(w)
wp = twist_log(exp_Result)
w

fprintf('Twist: Taylor test');
w = [rand(3,1); 0.05*rand(3,1)];
expm_Result = expm(se3_alg(w))
exp_Result  = twist_exp(w)
wp = twist_log(exp_Result)
w
sys = tf([0.1 8/15],[1 11/15 61/50 53/180 8/150]);
t = 0:0.7356:50;
u = sin(t) + 2*sin(2*t) + 3*sin(3*t) + 4*sin(4*t) + 5*sin(5*t) + 6*sin(6*t);
y = lsim(sys, u, t);

% Generating white noise with mean 0 and variance 1
noise = randn(size(y)); 
y_with_noise = y + noise;

subplot(2, 1, 1)
plot(t, u);
title("Input")

subplot(2, 1, 2)
plot(t, y_with_noise)
title("Output with White Noise")

% System identification
data = iddata(transpose(y_with_noise)', u', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Compare the estimated system with the true system
disp('True System:')
disp(sys);
disp('Estimated System:')
disp(estimated_sys);

% Calculate SSE
sse = sum((lsim(sys, u, t) - lsim(estimated_sys, u, t)).^2);
disp(['Sum of Squares Error (SSE): ', num2str(sse)]);

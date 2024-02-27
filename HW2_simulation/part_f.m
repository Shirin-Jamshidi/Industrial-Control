sys = tf([0.1 8/15],[1 11/15 61/50 53/180 8/150]);
t = 0:0.7356:50;

u6 = sin(t) + 2*sin(2*t) + 3*sin(3*t) + 4*sin(4*t) + 5*sin(5*t);
y6 = lsim(sys, u6, t);
subplot(6, 2, 1)
plot(t, y6);
title("y6");
subplot(6, 2, 2)
plot(t, u6);
title("u6");
% System identification
data = iddata(transpose(y6)', u6', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Calculate SSE
sse = sum((lsim(sys, u6, t) - lsim(estimated_sys, u6, t)).^2);
disp(['Sum of Squares Error (SSE): u6 ', num2str(sse)]);

u5 = sin(t) + 2*sin(2*t) + 3*sin(3*t) + 4*sin(4*t);
y5 = lsim(sys, u5, t);
subplot(6, 2, 3)
plot(t, y5);
title("y5");
subplot(6, 2, 4)
plot(t, u5);
title("u5");
% System identification
data = iddata(transpose(y5)', u5', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Calculate SSE
sse = sum((lsim(sys, u5, t) - lsim(estimated_sys, u5, t)).^2);
disp(['Sum of Squares Error (SSE): u5 ', num2str(sse)]);

u4 = sin(t) + 2*sin(2*t) + 3*sin(3*t);
y4 = lsim(sys, u4, t);
subplot(6, 2, 5)
plot(t, y4);
title("y4");
subplot(6, 2, 6)
plot(t, u4);
title("u4");
% System identification
data = iddata(transpose(y4)', u4', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Calculate SSE
sse = sum((lsim(sys, u4, t) - lsim(estimated_sys, u4, t)).^2);
disp(['Sum of Squares Error (SSE): u4 ', num2str(sse)]);

u3 = sin(t) + 2*sin(2*t);
y3 = lsim(sys, u3, t);
subplot(6, 2, 7)
plot(t, y3);
title("y3");
subplot(6, 2, 8)
plot(t, u3);
title("u4");
% System identification
data = iddata(transpose(y3)', u3', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Calculate SSE
sse = sum((lsim(sys, u3, t) - lsim(estimated_sys, u3, t)).^2);
disp(['Sum of Squares Error (SSE): u3 ', num2str(sse)]);

u2 = sin(t);
y2 = lsim(sys, u2, t);
subplot(6, 2, 9)
plot(t, y2);
title("y2");
subplot(6, 2, 10)
plot(t, u2);
title("u2");
% System identification
data = iddata(transpose(y2)', u2', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Calculate SSE
sse = sum((lsim(sys, u2, t) - lsim(estimated_sys, u2, t)).^2);
disp(['Sum of Squares Error (SSE): u2 ', num2str(sse)]);

u1 = t;
y1 = lsim(sys, u1, t);
subplot(6, 2, 11)
plot(t, y1);
title("y1");
subplot(6, 2, 12)
plot(t, u1);
title("u1");
% System identification
data = iddata(transpose(y1)', u1', 0.7356); % Create IDDATA object from input-output data
np = 4; 
nz = 1;

% Estimate transfer function parameters using least squares method
estimated_sys = tfest(data, np, nz, 'Ts', 0);

% Calculate SSE
sse = sum((lsim(sys, u1, t) - lsim(estimated_sys, u1, t)).^2);
disp(['Sum of Squares Error (SSE): u1 ', num2str(sse)]);
H = tf([0.1 8/15],[1 11/15 61/50 53/180 8/150]);
t = 0:0.7356:50;
u = sin(t) + 2*sin(2*t) + 3*sin(3*t) + 4*sin(4*t) + 5*sin(5*t) + 6*sin(6*t);
y = lsim(H, u, t); 
subplot(2, 1, 1)
plot(t, u);
title("output")
subplot(2, 1, 2)
plot(t, y)
title("input")
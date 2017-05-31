% error = [9.777; 0.707; 0.116; 0.029; 0.001]; % CN
error = [14.087; 5.71429205; 2.6000557; 1.24379; 0.61154831]; % Implicit
delta_t = [2.8097; 0.936567; 0.401386; 0.187313; 0.090635];
curve = fit(delta_t, error, 'power1')
figure(1)
plot(curve, delta_t, error)
figure(2)
loglog(delta_t, error)

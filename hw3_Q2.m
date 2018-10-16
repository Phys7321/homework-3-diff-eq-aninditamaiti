clear all;

dtheta0 = 0;
theta0_1 = 0.1;
theta0_2 = 0.2;
theta0_3 = 0.4;
theta0_4 = 0.8;
theta0_5 = 1;
w2 = 9;     % omega^2
g = 9.81;
R = g/w2;
m = 1;

% solutions of large oscillations

[T1, y1] = pendulum1(R,theta0_1, dtheta0, 1);
t1 = y1(:,1);
r1 = y1(:,2);  % theta
r1dot = y1(:,3); % d(theta)/dt
E1_total = 0.5*m*R^2*(r1dot.^2) + m*g*R.*(1 - cos(r1));

[T2, y2] = pendulum1(R,theta0_2, dtheta0, 1);
t2 = y2(:,1);
r2 = y2(:,2);  % theta
r2dot = y2(:,3); % d(theta)/dt
E2_total = 0.5*m*R^2*(r2dot.^2) + m*g*R.*(1 - cos(r2));

[T3, y3] = pendulum1(R,theta0_3, dtheta0, 1);
t3 = y3(:,1);
r3 = y3(:,2);  % theta
r3dot = y3(:,3); % d(theta)/dt
E3_total = 0.5*m*R^2*(r3dot.^2) + m*g*R.*(1 - cos(r3));

[T4, y4] = pendulum1(R,theta0_4, dtheta0, 1);
t4 = y4(:,1);
r4 = y4(:,2);  % theta
r4dot = y4(:,3); % d(theta)/dt
E4_total = 0.5*m*R^2*(r4dot.^2) + m*g*R.*(1 - cos(r4));

[T5, y5] = pendulum1(R,theta0_5, dtheta0, 1);
t5 = y5(:,1);
r5 = y5(:,2);  % theta
r5dot = y5(:,3); % d(theta)/dt
E5_total = 0.5*m*R^2*(r5dot.^2) + m*g*R.*(1 - cos(r5));

figure(1);
plot(t1,r1,'r-',t1,r1dot,'b-',t1,E1_total, 'g-');
legend('position for \theta_0=0.1','velocity for \theta_0=0.1','Total energy');
xlabel('t')
title('\theta_0 = 0.1')

figure(2);
plot(t2,r2,'r-',t2,r2dot,'b-',t2,E2_total, 'g-');
legend('position for \theta_0=0.2','velocity for \theta_0=0.2','Total energy');
xlabel('t')
title('\theta_0 = 0.2')

figure(3);
plot(t3,r3,'r-',t3,r3dot,'b-',t3,E3_total, 'g-');
legend('position for \theta_0=0.4','velocity for \theta_0=0.4','Total energy');
xlabel('t')
title('\theta_0 = 0.4')

figure(4);
plot(t4,r4,'r-',t4,r4dot,'b-',t4,E4_total, 'g-');
legend('position for \theta_0=0.8','velocity for \theta_0=0.8','Total energy');
xlabel('t')
title('\theta_0 = 0.8')

figure(5);
plot(t5,r5,'r-',t5,r5dot,'b-',t5,E5_total, 'g-');
legend('position for \theta_0=1','velocity for \theta_0=1','Total energy');
xlabel('t')
title('\theta_0 = 1')

sprintf('For each case both theta and d(theta)/dt are sinosoidal in nature and they differ by a phase difference of pi/2')
sprintf('For nonlinear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_1,T1,max(r1))
sprintf('For nonlinear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_2,T2,max(r2))
sprintf('For nonlinear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_3,T3,max(r3))
sprintf('For nonlinear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_4,T4,max(r4))
sprintf('For nonlinear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_5,T5,max(r5))
sprintf('Total energy is almost a constant for smaller amplitudes, as the amplitude keeps increasing, the variation in total energy keeps increasing too.')


T_array = [T1,T2,T3,T4,T5]; 
amp_array = [max(r1), max(r2), max(r3), max(r4), max(r5)];

% solutions by simple harmonic motion approximation

[T1p, y1p] = simple_pendulum(R,theta0_1, dtheta0, 1);
t1p = y1p(:,1);
r1p = y1p(:,2);  % theta
r1pdot = y1p(:,3); % d(theta)/dt

[T2p, y2p] = simple_pendulum(R,theta0_2, dtheta0, 1);
t2p = y2p(:,1);
r2p = y2p(:,2);  % theta
r2pdot = y2p(:,3); % d(theta)/dt

[T3p, y3p] = simple_pendulum(R,theta0_3, dtheta0, 1);
t3p = y3p(:,1);
r3p = y3p(:,2);  % theta
r3pdot = y3p(:,3); % d(theta)/dt

[T4p, y4p] = simple_pendulum(R,theta0_4, dtheta0, 1);
t4p = y4p(:,1);
r4p = y4p(:,2);  % theta
r4pdot = y4p(:,3); % d(theta)/dt

[T5p, y5p] = simple_pendulum(R,theta0_5, dtheta0, 1);
t5p = y5p(:,1);
r5p = y5p(:,2);  % theta
r5pdot = y5p(:,3); % d(theta)/dt

Tp_array = [T1p,T2p,T3p,T4p,T5p]; 
amp_p_array = [max(r1p), max(r2p), max(r3p), max(r4p), max(r5p)];

figure(6);
plot(amp_array,T_array, 'c-', amp_p_array, Tp_array, 'b-');
legend('nonlinear case', 'linear case');
title('T vs. \theta_{max}');
xlabel('\theta_{max}')
ylabel('T')

sprintf('For linear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_1,T1p,max(r1p))
sprintf('For linear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_2,T2p,max(r2p))
sprintf('For linear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_3,T3p,max(r3p))
sprintf('For linear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_4,T4p,max(r4p))
sprintf('For linear theta0 = %0.1f period is %0.3f and maximum amplitude is %0.3f',theta0_5,T5p,max(r5p))
sprintf('For nonlinear case as long as the amplitude is small, the time period remains constant. When amplitude increases beyond 0.8, time period increases too. This is because for large oscillations, the freqeuncy depends on the amplitude. The time period remains a constant for the linear case irrespective of amplitude. This is because the linear case is based on small amplitude approximation, hence the linear treatment is not correct for amplitude > 0.8. ')
sprintf('Time period in nonlinear case is larger. For nonlinear case, the frequency is determined by amplitude. So time period depends on amplitude too. Whereas in the linear case the frequency is a constant independent of the amplitude. Hence time period is a constant independent of the amplitude for linear case. ')

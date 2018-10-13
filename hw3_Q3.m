clear all;
theta_0 = 1; % initial displacement
dtheta_0 = 0;   % initial velocity
omega = 3;    % omega_2 = omega
g = 9.81;
R = g/omega^2;
m = 1;
gamma1 = 0.5;
gamma2 = 1;
gamma3 = 2;
gamma4 = 3;
gamma5 = 4;
gamma6 = 5;
gamma7 = 6;
gamma8 = 7;
gamma9 = 8;

[T1_graph, T1_theory , y1, te1] = pendulum2(R,gamma1 , theta_0, dtheta_0, 1);
t1 = y1(:,1);
r1 = y1(:,2);  % theta
r1dot = y1(:,3); % d(theta)/dt
omega1 = 2*pi/T1_graph;

figure(1);
plot(t1,r1,'k-',t1,r1dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 0.5')

[T2_graph, T2_theory , y2, te2] = pendulum2(R,gamma2 , theta_0, dtheta_0, 1);
t2 = y2(:,1);
r2 = y2(:,2);  % theta
r2dot = y2(:,3); % d(theta)/dt
omega2 = 2*pi/T2_graph;

figure(2);
plot(t2,r2,'k-',t2,r2dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 1')

[T3_graph, T3_theory, y3, te3] = pendulum2(R,gamma3 , theta_0, dtheta_0, 1);
t3 = y3(:,1);
r3 = y3(:,2);  % theta
r3dot = y3(:,3); % d(theta)/dt
omega3 = 2*pi/T3_graph;

figure(3);
plot(t3,r3,'k-',t3,r3dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 2')

[T4_graph, T4_theory, y4, te4] = pendulum2(R,gamma4 , theta_0, dtheta_0, 1);
t4 = y4(:,1);
r4 = y4(:,2);  % theta
r4dot = y4(:,3); % d(theta)/dt
omega4 = 2*pi/T4_theory;

figure(4);
plot(t4,r4,'k-',t4,r4dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 3')

[T1_undamped , y_undamped] = pendulum1(R, theta_0, dtheta_0, 1);
t_undamped = y_undamped(:,1);
r_undamped = y_undamped(:,2);  % theta
r_undampeddot = y_undamped(:,3); % d(theta)/dt
omega_undamped = 2*pi/T1_undamped;

sprintf('For gamma = %0.1f period is %0.3f and angular frequency is %0.3f',gamma1, T1_graph, omega1)
sprintf('For gamma = %0.0f period is %0.3f and angular frequency is %0.3f',gamma2, T2_graph, omega2)
sprintf('For gamma = %0.0f period is %0.3f and angular frequency is %0.3f',gamma3, T3_graph, omega3)
sprintf('For gamma = %0.0f period is %0.3f and angular frequency is %0.3f',gamma4, T4_theory, omega4)
sprintf('Gamma = 3 is critical damping hence T = 0')
sprintf('For undamped case, period is %0.3f and angular frequency is %0.3f', T1_undamped, omega_undamped)
sprintf('Theoretical angular frequency keeps decreasing with increasing damping until critical damping is reached.')
sprintf('Can be seen from theoretical T = %0.3f for gamma = %0.0f', T1_theory, gamma1)
sprintf('Can be seen from theoretical T = %0.3f for gamma = %0.0f', T2_theory, gamma2)
sprintf('Can be seen from theoretical T = %0.3f for gamma = %0.0f', T3_theory, gamma3)
sprintf('As theoretically T increases, frequency and angular frequency decreases')
sprintf('Period of gamma = 0.5 is lesser than undamped case')

[T5_graph, T5_theory, y5, te5] = pendulum2(R,gamma5 , theta_0, dtheta_0, 1);
t5 = y5(:,1);
r5 = y5(:,2);  % theta
r5dot = y5(:,3); % d(theta)/dt

figure(5);
plot(t5,r5,'k-',t5,r5dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 4')

[T6_graph, T6_theory, y6, te6] = pendulum2(R,gamma6 , theta_0, dtheta_0, 1);
t6 = y6(:,1);
r6 = y6(:,2);  % theta
r6dot = y6(:,3); % d(theta)/dt

figure(6);
plot(t6,r6,'k-',t6,r6dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 5')

[T7_graph, T7_theory, y7, te7] = pendulum2(R,gamma7 , theta_0, dtheta_0, 1);
t7 = y7(:,1);
r7 = y7(:,2);  % theta
r7dot = y7(:,3); % d(theta)/dt

figure(7);
plot(t7,r7,'k-',t7,r7dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 6')

[T8_graph, T8_theory, y8, te8] = pendulum2(R,gamma8 , theta_0, dtheta_0, 1);
t8 = y8(:,1);
r8 = y8(:,2);  % theta
r8dot = y8(:,3); % d(theta)/dt

figure(8);
plot(t8,r8,'k-',t8,r8dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 7')

[T9_graph, T9_theory, y9, te9] = pendulum2(R,gamma9 , theta_0, dtheta_0, 1);
t9 = y9(:,1);
r9 = y9(:,2);  % theta
r9dot = y9(:,3); % d(theta)/dt

figure(9);
plot(t9,r9,'k-',t9,r9dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 8')

sprintf('For gamma > omega the motion is not oscillatory, it''s exponentially decreasing motion')
sprintf('As gamma increases, it takes more time to reach equilibrium')
sprintf('Gamma = %0.0f decays to equilibrium at t = %0.3f', gamma4, te4)
sprintf('Gamma = %0.0f decays to equilibrium at t = %0.3f', gamma5, te5)
sprintf('Gamma = %0.0f decays to equilibrium at t = %0.3f', gamma6, te6)
sprintf('Gamma = %0.0f decays to equilibrium at t = %0.3f', gamma7, te7)
sprintf('Gamma = %0.0f decays to equilibrium at t = %0.3f', gamma8, te8)
sprintf('Gamma = %0.0f decays to equilibrium at t = %0.3f', gamma9, te9)

figure(10)
plot(r9,r9dot,'b-', r7, r7dot, 'r-', r5, r5dot, 'g-', r3, r3dot, 'k-', r1, r1dot, 'c-')
legend('\gamma=8','\gamma=6','\gamma=4','\gamma=2','\gamma=0.5')
title('Phase space diagrams')

sprintf('Qualitative feature of each phase space diagram depends on gamma. For gamma>omega corresponding to overdamped motion, the phase space diagrams are similar. The one for critical damping is distinct. Those for under-damped motion are similar.')

clear all;
theta1_0 = 1; % initial displacement
dtheta1_0 = 0;   % initial velocity
omega0 = 3;    % omega_2 = omega
g = 9.81;
R = g/omega0^2;
m = 1;
gamma1 = 0.5;
omegap1 = 2;
A0 = 1;

[T1 , y1, te1] = driven_damped(R,gamma1 , theta1_0, dtheta1_0, 1, omegap1, A0);
t1 = y1(:,1);
r1 = y1(:,2);  % theta
r1dot = y1(:,3); % d(theta)/dt
omega1 = 2*pi/T1;

[T1_graph, T1_theory , y1_u, te1_u] = pendulum2(R,gamma1 , theta1_0, dtheta1_0, 1);
t1_u = y1_u(:,1);
r1_u = y1_u(:,2);  % theta
r1dot_u = y1_u(:,3); % d(theta)/dt
omega1_u = 2*pi/T1_graph;

figure(1);
plot(t1,r1,'k-',t1,r1dot,'b-', t1_u, r1_u, 'r-', t1_u, r1dot_u, 'g-')
legend('Position of driven','Velocity of driven','Position of undriven','Velocity of undriven')
title('Position and velocity vs time for driven & undriven case with , \theta_0 = 1, d\theta_0/dt = 0')

sprintf('The driven case oscillates with constant frequency after a while.')
sprintf('Period is %0.3f and frequency is %0.3f after several oscillations', T1, omega1)

theta2_0 = 0; % initial displacement
dtheta2_0 = 1;   % initial velocity

[T2 , y2, te2] = driven_damped(R,gamma1 , theta2_0, dtheta2_0, 1, omegap1, A0);
t2 = y2(:,1);
r2 = y2(:,2);  % theta
r2dot = y2(:,3); % d(theta)/dt
omega2 = 2*pi/T2;

figure(2);
plot(t2,r2,'k-',t2,r2dot,'b-')
legend('Position of driven','Velocity of driven')
title('Position and velocity vs time for \omega=2, \theta_0 = 0, d\theta_0/dt = 1')

sprintf('The driven case oscillates with constant frequency after a while.')
sprintf('Period is %0.3f and frequency is %0.3f after several oscillations', T2, omega2)

sprintf('It can be seen that the first part of the oscilllations is transient as it depends on initial conditions and decays with time. After the transient part a steady oscillation can be observed which is independent of the initial conditions. This is due to the driving force and the amplitude and frequency of the driving force are what determines the stead part')
sprintf('theta(t) approiaches a limiting behavior independent of initial conditions, this behavior is determined by A0 and omega')

omegap2 = 1;
omegap3 = 4;

[T3 , y3, te3] = driven_damped(R,gamma1 , theta1_0, dtheta1_0, 1, omegap2, A0);
t3 = y3(:,1);
r3 = y3(:,2);  % theta
r3dot = y3(:,3); % d(theta)/dt
omega3 = 2*pi/T3;

figure(3);
plot(t3,r3,'k-',t3,r3dot,'b-')
legend('Position of driven','Velocity of driven')
title('Position and velocity vs time for \omega=1, \theta_0 = 1, d\theta_0/dt = 0')
sprintf('Period is %0.3f and frequency is %0.3f after several oscillations', T3, omega3)


[T4 , y4, te4] = driven_damped(R,gamma1 , theta2_0, dtheta2_0, 1, omegap2, A0);
t4 = y4(:,1);
r4 = y4(:,2);  % theta
r4dot = y4(:,3); % d(theta)/dt
omega4 = 2*pi/T4;

figure(4);
plot(t4,r4,'k-',t4,r4dot,'b-')
legend('Position of driven','Velocity of driven')
title('Position and velocity vs time for \omega=1, \theta_0 = 0, d\theta_0/dt = 1')
sprintf('Period is %0.3f and frequency is %0.3f after several oscillations', T4, omega4)

[T5 , y5, te5] = driven_damped(R,gamma1 , theta1_0, dtheta1_0, 1, omegap3, A0);
t5 = y5(:,1);
r5 = y5(:,2);  % theta
r5dot = y5(:,3); % d(theta)/dt
omega5 = 2*pi/T5;

figure(5);
plot(t5,r5,'k-',t5,r5dot,'b-')
legend('Position of driven','Velocity of driven')
title('Position and velocity vs time for \omega=4, \theta_0 = 1, d\theta_0/dt = 0')
sprintf('Period is %0.3f and frequency is %0.3f after several oscillations', T5, omega5)


[T6 , y6, te6] = driven_damped(R,gamma1 , theta2_0, dtheta2_0, 1, omegap3, A0);
t6 = y6(:,1);
r6 = y6(:,2);  % theta
r6dot = y6(:,3); % d(theta)/dt
omega6 = 2*pi/T6;

figure(6);
plot(t6,r6,'k-',t6,r6dot,'b-')
legend('Position of driven','Velocity of driven')
title('Position and velocity vs time for \omega=4, \theta_0 = 0, d\theta_0/dt = 1')
sprintf('Period is %0.3f and frequency is %0.3f after several oscillations', T6, omega6)


sprintf('As can be seen from the plots, the frequency the steady part doubles for omega=4')
sprintf('The frequency the steady part halvs for omega=1. Omega = frequency of driving force')
sprintf('The amplitude of steady part doesn''t depend on initial conditions')
sprintf('Amplitude of steady part is determined by A0 and omega of driving force')

wp = [0, 1.0, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4];
delta1 = [];
for i=1:length(wp)
    omega_prime = wp(i);
    [Tp , yp, tep] = driven_damped(R,gamma1 , theta1_0, dtheta1_0, 1, omega_prime, A0);
    tp = yp(:,1);
    rp = yp(:,2);  % theta
    rpdot = yp(:,3); % d(theta)/dt
    
    A = A0/((omega0^2 - omega_prime^2)^2 + 4*(gamma1*omega_prime)^2)^0.5;
    y1_steady = A.*cos(omega_prime.* tp);  % steady state solution
    t_length = length(tp);
    t1_length = int16(t_length/2);
    delta =rp(end-t1_length:end)-y1_steady(end-t1_length:end);   
    delta = mean(delta);    % average differenc between applied forced motion and steady state motion
    delta1 = [delta1, delta];
    y1_steady = y1_steady + delta;
    figure(6+i)
    plot(tp,y1_steady, 'k',tp,rp,'g');
    xlabel('t')
    legend('Steady solution','Driven solution')
    title('Phase shift in steady and driven solution for \gamma=0.5')
    
end

gamma2 = 1.5;
delta2 = [];
for i=1:length(wp)
    omega_prime = wp(i);
    [Tp , yp, tep] = driven_damped(R, gamma2 , theta1_0, dtheta1_0, 1, omega_prime, A0);
    tp = yp(:,1);
    rp = yp(:,2);  % theta
    rpdot = yp(:,3); % d(theta)/dt
    
    A = A0/((omega0^2 - omega_prime^2)^2 + 4*(gamma2*omega_prime)^2)^0.5;
    y1_steady = A.*cos(omega_prime.* tp);  % steady state solution
    t_length = length(tp);
    t1_length = int16(t_length/2);
    delta =rp(end-t1_length:end)-y1_steady(end-t1_length:end);   
    delta = mean(delta);%/pi;    % average differenc between applied forced motion and steady state motion
    delta2 = [delta2, delta];
    y1_steady = y1_steady + delta; %A.*cos(omega_prime.* tp + phase);
    figure(16+i)
    plot(tp,y1_steady, 'k',tp,rp,'g');
    xlabel('t')
    legend('Steady solution','Driven solution')
    title('Phase shift in steady and driven solution for \gamma=1.5')
    
end

figure(27);
plot(wp, delta1, 'c--', wp,delta2, 'm-');
legend('\gamma = 0.5','\gamma=1.5');
xlabel('\omega');
ylabel('\delta');
title('\delta(\omega)')

sprintf('For smaller omega, the difference between delta values are smaller. As omega grows, the difference between delta of two gammas increase.')



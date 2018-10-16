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

KE1 = 0.5*m*R^2.*r1dot.*r1dot;
PE1 = m*g*R*(1-cos(r1));
E1total =KE1 + PE1;
KE1_avg =[];
PE1_avg =[];
E1_avg =[];
ind1 = find(r1dot.*circshift(r1dot , [-1 0]) <= 0);
for i =1:(length(ind1)-2) / 2 
    n1 =2 .* i -1;
    KE_element = 0.5 .* R.^2 .* r1dot(ind1(n1): ind1(n1+2)).^2;  
    PE_element = g .* R .* (1-cos(r1(ind1(n1): ind1(n1)))) ; 
    E_element=KE_element + PE_element; 
    KE1_avg =[KE1_avg, mean(KE_element)];
    PE1_avg =[PE1_avg, mean(PE_element)];
    E1_avg =[E1_avg, mean(E_element)];
    
end    
deltan1 =1:(length(ind1)-2) / 2 ;

figure(1);
plot(t1,r1,'k-',t1,r1dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 0.5')

figure(2);
plot(deltan1,KE1_avg,deltan1,PE1_avg,deltan1,E1_avg);
legend('Average potential energy', 'Avergae kinetic energy', 'Average total energy');
xlabel('Number of cycles n');
ylabel('Energies');
title('Average energies as a function of number of cycles for \gamma = 0.5')


[T2_graph, T2_theory , y2, te2] = pendulum2(R,gamma2 , theta_0, dtheta_0, 1);
t2 = y2(:,1);
r2 = y2(:,2);  % theta
r2dot = y2(:,3); % d(theta)/dt
omega2 = 2*pi/T2_graph;

KE2 = 0.5*m*R^2.*r2dot.*r2dot;
PE2 = m*g*R*(1-cos(r2));
E2total =KE2 + PE2;
KE2_avg =[];
PE2_avg =[];
E2_avg =[];
ind2 = find(r2dot.*circshift(r2dot , [-1 0]) <= 0);
for i =1:(length(ind2)-2) / 2 
    n1 =2 .* i -1;
    KE_element = 0.5 .* R.^2 .* r2dot(ind2(n1): ind2(n1+2)).^2;  
    PE_element = g .* R .* (1-cos(r2(ind2(n1): ind2(n1)))) ; 
    E_element=KE_element + PE_element; 
    KE2_avg =[KE2_avg, mean(KE_element)];
    PE2_avg =[PE2_avg, mean(PE_element)];
    E2_avg =[E2_avg, mean(E_element)];
    
end    
deltan2 =1:(length(ind2)-2) / 2 ;


figure(3);
plot(t2,r2,'k-',t2,r2dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 1')

figure(4);
plot(deltan2,KE2_avg,deltan2,PE2_avg,deltan2,E2_avg);
legend('Average potential energy', 'Avergae kinetic energy', 'Average total energy');
xlabel('Number of cycles n');
ylabel('Energies');
title('Average energies as a function of number of cycles for \gamma = 1')

[T3_graph, T3_theory, y3, te3] = pendulum2(R,gamma3 , theta_0, dtheta_0, 1);
t3 = y3(:,1);
r3 = y3(:,2);  % theta
r3dot = y3(:,3); % d(theta)/dt
omega3 = 2*pi/T3_graph;

KE3 = 0.5*m*R^2.*r3dot.*r3dot;
PE3 = m*g*R*(1-cos(r3));
E3total =KE3 + PE3;
KE3_avg =[];
PE3_avg =[];
E3_avg =[];
ind3 = find(r3dot.*circshift(r3dot , [-1 0]) <= 0);
for i =1:(length(ind3)-2) / 2 
    n1 =2 .* i -1;
    KE_element = 0.5 .* R.^2 .* r3dot(ind3(n1): ind3(n1+2)).^2;  
    PE_element = g .* R .* (1-cos(r3(ind3(n1): ind3(n1)))) ; 
    E_element=KE_element + PE_element; 
    KE3_avg =[KE3_avg, mean(KE_element)];
    PE3_avg =[PE3_avg, mean(PE_element)];
    E3_avg =[E3_avg, mean(E_element)];
    
end    
deltan3 =1:(length(ind3)-2) / 2 ;

figure(5);
plot(t3,r3,'k-',t3,r3dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 2')

figure(6);
plot(deltan3,KE3_avg,deltan3,PE3_avg,deltan3,E3_avg);
legend('Average potential energy', 'Avergae kinetic energy', 'Average total energy');
xlabel('Number of cycles n');
ylabel('Energies');
title('Average energies as a function of number of cycles for \gamma = 2')

[T4_graph, T4_theory, y4, te4] = pendulum2(R,gamma4 , theta_0, dtheta_0, 1);
t4 = y4(:,1);
r4 = y4(:,2);  % theta
r4dot = y4(:,3); % d(theta)/dt
omega4 = 2*pi/T4_graph;

KE4 = 0.5*m*R^2.*r4dot.*r4dot;
PE4 = m*g*R*(1-cos(r4));
E4total =KE4 + PE4;
KE4_avg =[];
PE4_avg =[];
E4_avg =[];
ind4 = find(r4dot.*circshift(r4dot , [-1 0]) <= 0);
for i =1:(length(ind4)-2) / 2 
    n1 =2 .* i -1;
    KE_element = 0.5 .* R.^2 .* r4dot(ind4(n1): ind4(n1+2)).^2;  
    PE_element = g .* R .* (1-cos(r4(ind4(n1): ind4(n1)))) ; 
    E_element=KE_element + PE_element; 
    KE4_avg =[KE4_avg, mean(KE_element)];
    PE4_avg =[PE4_avg, mean(PE_element)];
    E4_avg =[E4_avg, mean(E_element)];
    
end    
deltan4 =1:(length(ind4)-2) / 2 ;

figure(7);
plot(t4,r4,'k-',t4,r4dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 3')

figure(8);
plot(t4,KE4,'c-',t4,PE4,'m-',t4,E4total,'y-', deltan4, KE4_avg, 'k-', deltan4, PE4_avg, 'k-', deltan4, E4_avg, 'k-');
legend('Potential energy', 'Kinetic energy', 'Total energy');
xlabel('Time');
ylabel('Energies');
title('Energies as a function of time for \gamma = 3')

[T1_undamped , y_undamped] = pendulum1(R, theta_0, dtheta_0, 1);
t_undamped = y_undamped(:,1);
r_undamped = y_undamped(:,2);  % theta
r_undampeddot = y_undamped(:,3); % d(theta)/dt
omega_undamped = 2*pi/T1_undamped;
gamma_matrix = [0, gamma1, gamma2, gamma3, gamma4];
omega_matrix = [omega_undamped, omega1, omega2, omega3, omega4];
T_matrix = [T1_undamped, T1_graph, T2_graph, T3_graph, T4_graph];

figure(9);
plot(gamma_matrix, omega_matrix, 'c-', gamma_matrix, T_matrix, 'm-');
legend('\omega', 'T');
xlabel('\gamma');
title('\omega and T as function of \gamma')


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

figure(10);
plot(t5,r5,'k-',t5,r5dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 4')

[T6_graph, T6_theory, y6, te6] = pendulum2(R,gamma6 , theta_0, dtheta_0, 1);
t6 = y6(:,1);
r6 = y6(:,2);  % theta
r6dot = y6(:,3); % d(theta)/dt

figure(11);
plot(t6,r6,'k-',t6,r6dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 5')

[T7_graph, T7_theory, y7, te7] = pendulum2(R,gamma7 , theta_0, dtheta_0, 1);
t7 = y7(:,1);
r7 = y7(:,2);  % theta
r7dot = y7(:,3); % d(theta)/dt

figure(12);
plot(t7,r7,'k-',t7,r7dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 6')

[T8_graph, T8_theory, y8, te8] = pendulum2(R,gamma8 , theta_0, dtheta_0, 1);
t8 = y8(:,1);
r8 = y8(:,2);  % theta
r8dot = y8(:,3); % d(theta)/dt

figure(13);
plot(t8,r8,'k-',t8,r8dot,'b-')
legend('Position','Velocity')
title('Position and velocity vs time for \gamma = 7')

[T9_graph, T9_theory, y9, te9] = pendulum2(R,gamma9 , theta_0, dtheta_0, 1);
t9 = y9(:,1);
r9 = y9(:,2);  % theta
r9dot = y9(:,3); % d(theta)/dt

figure(14);
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

figure(15)
plot(r9,r9dot,'b-', r7, r7dot, 'r-', r5, r5dot, 'g-', r3, r3dot, 'k-', r1, r1dot, 'c-')
legend('\gamma=8','\gamma=6','\gamma=4','\gamma=2','\gamma=0.5')
title('Phase space diagrams')

sprintf('Qualitative feature of each phase space diagram depends on gamma. For gamma>omega corresponding to overdamped motion, the phase space diagrams are similar. The one for critical damping is distinct. Those for under-damped motion are similar.')

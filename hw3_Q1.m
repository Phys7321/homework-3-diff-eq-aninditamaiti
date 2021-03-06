clear all;

theta_0 = 0.25; % initial displacement
dtheta_0 = 0;   % initial velocity
omega_2 = 3;    % omega_2 = omega
m = 1;
g = 9.81;
R = g/omega_2^2;

[T, y] = simple_pendulum(R,theta_0, dtheta_0, 1);
t = y(:,1);
r = y(:,2);  % theta
rdot = y(:,3); % d(theta)/dt

figure(1);
plot(t,r,'k-',t,rdot,'b-')
legend('Position','Velocity')
xlabel('t')
title('Position and velocity vs time')


E_K =  0.5*m*R^2*(rdot.^2);
E_P = 0.5*m*R^2*( omega_2^2*r.^2);
E_total = E_K + E_P;

figure(2);
plot(t,E_total,'c-', t, E_K, 'm-', t, E_P, 'k-')    % plotting total energy vs time
legend('Total energy', 'Kinetic energy', 'Potential energy')
xlabel('t')
title('Energies as a function of time')

period = 2*pi/omega_2;
index = [];
t_delta_n =[];
delta_n = [];
for i= 1: size(t)-1
    for j=1:10
        if (abs(t(i) - j*period)<0.0225)
            index = [index;i];
            t_delta_n = [t_delta_n;t(i)];
            delta_i = (E_total(i)-E_total(1))/E_total(1);
            delta_n = [delta_n;delta_i];
        end
    end
end

delta_n;
figure(3);
plot(t_delta_n,delta_n,'-')    % plotting relative chaneg in total energy for each cycle
legend('\delta E_n')
ylabel('\delta E_n');
xlabel('n')
title('relative change in total energy during each cycle')
sprintf("Relative change in total energy not unifrom over each cycle")


delta = [];
for i=1:size(E_total)
    delta_i = (E_total(i)-E_total(1))/E_total(1);
    delta = [delta;delta_i];
end
figure(4);
plot(t,delta,'-')    % plotting relative chaneg in total energy vs time
legend('\delta E')
xlabel('t')
ylabel('\delta E')
title('relative change in total energy vs time')
sprintf("Relative change in total energy not unifrom over time")
    
sum_KE = 0;
sum_PE = 0;


for i=index(1):index(2);
    sum_PE = sum_PE + E_P(i); 
    sum_KE = sum_KE + E_K(i); 
end
avg_PE = sum_PE/period;
avg_KE = sum_KE/period;
t_avg = t(index(1):index(2));

a_PE=[]; a_KE=[];
for i=index(1):index(2)
    a_PE = [a_PE;avg_PE];
    a_KE = [a_KE;avg_KE];
end



KE_avg =[];
PE_avg =[];
E_avg =[];
ind = find(rdot.*circshift(rdot , [-1 0]) <= 0);
for i =1:(length(ind)-2) / 2 
    n1 =2 .* i -1;
    KE_element = 0.5 .* R.^2 .* rdot(ind(n1): ind(n1+2)).^2;  
    PE_element = 0.5 .* m .* R^2.*( omega_2^2*r(ind(n1): ind(n1+2)).^2);%g .* R .* (1-cos(r(ind(n1): ind(n1 +2)))) ; 
    E_element=KE_element + PE_element; 
    KE_avg =[KE_avg, mean(KE_element)];
    PE_avg =[PE_avg, mean(PE_element)];
    E_avg =[E_avg, mean(E_element)];
    
end    
n  =1:(length(ind)-2) / 2 ;




figure(5);
plot(n,PE_avg,'b-',n,KE_avg,'k-')    % plotting average potential & kientic energy vs time
legend('Avergae Potential Energy','Average Kinetic Energy')
title('Average potential and kinetic energy over one cycle')
xlabel('n')
ylabel('Average energies')

figure(6);
plot(r, rdot,'k-')    % phase space
legend('\omega^2 = 9')
title('Phase space diagram')
xlabel('\theta')
ylabel('d\theta / dt')


theta_01 = 0; % initial displacement
dtheta_01 = 1; % initial velocity

[T1, y1] = simple_pendulum(R,theta_01, dtheta_01, 1);
r1 = y1(:,2);  % theta
r1dot = y1(:,3); % d(theta)/dt

figure(7);
plot(r1, r1dot,'k-')    % phase space
legend('\omega^2 = 9')
title('Phase space diagram')
xlabel('\theta')
ylabel('d\theta / dt')

theta_02 = 1; % initial displacement
dtheta_02 = 1; % initial velocity

[T2, y2] = simple_pendulum(R,theta_02, dtheta_02, 1);
r2 = y2(:,2);  % theta
r2dot = y2(:,3); % d(theta)/dt

figure(8);
plot(r2, r2dot,'k-')    % phase space
legend('omega^2 = 9')
title('Phase space diagram')
xlabel('\theta')
ylabel('d\theta / dt')

sprintf("Depending on initial displacement and initial velocity all paths are different.")
sprintf("Shape of all paths are same, i.e. ellipse.") 
sprintf("Each path is characterized by unique initial conditions.")
sprintf("Motion of a representative point (theta, d(theta)/dt is always counter-clockwise.")
sprintf("This can be understood by plotting first 10, first 20, first 30 and so on representative points for each case.")
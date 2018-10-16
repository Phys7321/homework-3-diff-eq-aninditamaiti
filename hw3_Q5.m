clear all;
theta1_0 = 0; % initial displacement
dtheta1_0 = 0;   % initial velocity
omega0 = 3;    % omega_2 = omega
g = 9.81;
R = g/omega0^2;
m = 1;
gamma = [0.1, 0.5,1,2];
omegap = [0, 1.0, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4];
A0 = 1;
Amax_matrix = [];
Q_inverse_matrix = [];

for j=1:length(gamma)
    gamma1 = gamma(j);
    Amp = [];

    for i=1: length(omegap)
        
        omega = omegap(i);
        A = A0/((omega0^2 - omega^2)^2 + 4*(gamma1*omega)^2)^0.5;
        Amp = [Amp,A];
    end

    figure(j)
    plot(omegap, Amp, 'b')
    xlabel('\omega')
    ylabel('Amplitude(\omega)')
    title('Amplitude vs \omega for different \gamma')

    [A_max, index]= max(Amp);
    w_max = omegap(index);
    sprintf('Maximum amplitude %0.0f at w=%0.1f', A_max, w_max)
    sprintf('W_max is same as natural frequency w0')

    A_half = A_max/sqrt(2)*ones(1,length(Amp));
    diff = abs(Amp - A_half);
    diff = sort(diff);
    min2 = [diff(1), diff(2)];
    index3 = [];
    for i=1:length(Amp)
        if abs(abs(Amp(i)-A_half)-min2(1))<eps | abs(abs(Amp(i)-A_half)-min2(2))<eps
            index3 = [index3, i];
        end
    end

    delta_w = omegap(index3(2)) - omegap(index3(1));
    Q_inverse = delta_w/w_max;
    Amax_matrix = [Amax_matrix, A_max];
    Q_inverse_matrix = [Q_inverse_matrix, Q_inverse];


end

figure(length(gamma)+1)
plot(gamma, Amax_matrix, 'b-', gamma, Q_inverse_matrix, 'r-')
legend('A_{max}','\Delta\omega/\omega_{max}')
xlabel('\gamma')

sprintf('A_max decreases with increasing gamma')
sprintf('1/Q reaches a peak at gamma=1, i.e. quality factor is lowest at this point. A local maxima.')
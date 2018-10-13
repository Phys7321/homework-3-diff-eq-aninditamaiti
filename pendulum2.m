function [period,T, sol,te] = pendulum2(R, gamma ,theta0,thetad0,grph) 
% Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0
    error('Must input length and initial conditions')
end
if nargin==1
   theta0 = pi/2;
   thetad0=0;
   grph=0;
end
if nargin==2
    thetad0 = 0;
    grph=1;
end
if nargin==3
    grph=1;
end
g=9.81;
omega2 = sqrt(g/R);
omega = sqrt(omega2^2 - gamma^2);


if omega>0
    T= 2*pi/omega;
    % number of oscillations to graph
    N = 10;
    tspan = [0 N*T];
    opts = odeset('refine',6);
else
    T = 0;
    tspan = [0 40];
    opts = odeset('Events',@stopevents);
end

%opts = odeset('events',@events,'refine',6); %Here for future event finder

r0 = [theta0 thetad0];
[t,w, te] = ode45(@proj,tspan,r0,opts,g,R, gamma);
sol = [t,w];
ind= find(w(:,2).*circshift(w(:,2), [-1 0]) <= 0);
ind = chop(ind,4);
period= 2*mean(diff(t(ind)));



end
%-------------------------------------------
%
function rdot = proj(t,r,g,R, gamma)
    rdot = [r(2); (-2*gamma*r(2)- g*sin(r(1))/R)];
end

function [position,isterminal,direction] = stopevents(t,r,g,R,gamma)
    position = r(1)-0.000099; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end
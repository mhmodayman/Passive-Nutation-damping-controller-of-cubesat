%close all;
clear, clc;

%--------------------------------------------------------------------------
sc = 1; % = 1 for mercury and 2 for oil
%--------------------------------------------------------------------------

global g R dR H D gama J J3 rho nu K Kd sigma m cd I1 I2 I3 I4

g = 9.8;                    % m/s^2
R = 16.71e-3;               % m
dR = 1.4e-3;                % m
H = 30.39e-3;               % m
D = 0.78e-3;                % m
gama = 1.33;                % rad
J = 4000;                   % kg.m^2
J3 = 6000;                  % kg.m^2
rho = [13600 912];          % kg/m^3    mercury -- oil
nu = [1.17e-7 0.00042];     % m^2/s    mercury -- oil

x1 = 1;         % e0
x2 = 0;         % e1
x3 = 0;         % e2
x4 = 0;         % e3
x5 = 100;       % omega_x
x6 = 0;         % omega_y
x7 = 250*pi;    % omega_z
x8 = 0;         % beta_ dot
x9 = 0;         % beta

K = sin(gama)/gama;
Kd = (sin(gama/2))*(2/gama);
sigma = J3/J;
m = gama*R*D*rho(sc)*dR;
cd = 0.0665*rho(sc)*(D+dR)*R*gama*(((sigma-1)*R*x7)^0.75)*((((D+dR)*nu(sc))/(D*dR))^0.25) +...
     0.00593*rho(sc)*(((D+dR)*D*dR)^0.5)*(R^1.5)*gama*(sigma-1)*x7;


I1 = m*(H*H + R*R*(1-K)*.5);
I2 = m*(H*H + R*R*(1+K)*.5);
I3 = m*R*R;
I4 = m*R*H*Kd;

x0 = [x1 x2 x3 x4 x5 x6 x7 x8 x9]';   % initial conditions

tspan = [0 1];

[t,x] = ode45(@(t,x) odefc(x), tspan, x0);

qq = x(:,1:4);

t = 500*t;

ang = ones(3,1)*0.1;

for ii = 2 : length(x)
    ang(:,ii) = angls(qq(ii-1,:), ang(:,ii-1));
end
ang = ang';

th = rad2deg(ang);
plot(t, th(:,3)), 
if sc == 1
    title('mercury')
else
    title('oil')
end
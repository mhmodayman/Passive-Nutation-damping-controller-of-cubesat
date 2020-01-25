function gg = odefc(x)

global g R dR H D gama J J3 rho nu K Kd sigma m cd I1 I2 I3 I4

f1 = (J-I3-J3+I1*sin(x(9))*sin(x(9))+I2*cos(x(9))*cos(x(9)))*x(6)*x(7) +...
     I4*cos(x(9))*x(5)*x(6) + (I1-I2)*sin(2*x(9))*x(5)*(x(7)/2) + I4*sin(x(9))*x(6)*x(6) -...
     I4*sin(x(9))*x(7)*x(7) + ((I1-I2)*sin(2*x(9)))*x(5)*x(8) -...
     2*I4*sin(x(9))*x(7)*x(8) + ((I2-I1)*cos(2*x(9))-I3)*x(6)*x(8) - I4*sin(x(9))*x(8)*x(8) +...
     2*m*g*H*(x(1)*x(2)+x(3)*x(4)) + m*g*Kd*R*(-x(1)^2+x(2)^2+x(3)^2-x(4)^2)*sin(x(9));

f2 = (I3-J+J3-I2*sin(x(9))*sin(x(9))-I1*cos(x(9))*cos(x(9)))*x(5)*x(7) - I4*sin(x(9))*x(5)*x(6) +...
     (I2-I1)*sin(2*x(9))*x(6)*(x(7)/2) - I4*cos(x(9))*x(5)*x(5) +...
     I4*cos(x(9))*x(7)*x(7) + I4*cos(x(9))*x(8)*x(8) + (I3+(I2-I1)*cos(2*x(9)))*x(5)*x(8) + 2*I4*cos(x(9))*x(7)*x(8) +...
     (I2-I1)*sin(2*x(9))*x(6)*x(8) + 2*m*g*H*(x(1)*x(3)-x(2)*x(4)) + m*g*Kd*R*(x(1)^2-x(2)^2-x(3)^2+x(4)^2)*cos(x(9));

f3 = (I1-I2)*cos(2*x(9))*x(5)*x(6) + I4*sin(x(9))*x(5)*x(7) - I4*cos(x(9))*x(6)*x(7) +...
     (I1-I2)*sin(2*x(9))*((x(6)*x(6)-x(5)*x(5))/2) - 2*m*g*Kd*R*(x(1)*x(2)+x(3)*x(4))*cos(x(9)) +...
     2*m*g*Kd*R*(x(2)*x(4)-x(1)*x(3))*sin(x(9));

f4 = (I1-I2)*cos(2*x(9))*x(5)*x(6) - I4*cos(x(9))*x(6)*x(7) + I4*sin(x(9))*x(5)*x(7) +...
     ((I1-I2)*sin(2*x(9))*(x(6)*x(6)-x(5)*x(5))/2) - 2*m*g*Kd*R*(x(1)*x(2)+x(3)*x(4))*cos(x(9)) +...
     2*m*g*Kd*R*(x(2)*x(4)-x(1)*x(3))*sin(x(9))-cd*R*R*x(8);


M = [J+I1*cos(x(9))*cos(x(9))+I2*sin(x(9))*sin(x(9))     (I1-I2)*sin(2*x(9))/2                               -I4*cos(x(9))     -I4*cos(x(9))     0
     (I1-I2)*sin(2*x(9))/2                               J+I2*cos(x(9))*cos(x(9))+I1*sin(x(9))*sin(x(9))     -I4*sin(x(9))     -I4*sin(x(9))     0
     -I4*cos(x(9))                                       -I4*sin(x(9))                                       I3+J3             I3                0
     -I4*cos(x(9))                                       -I4*sin(x(9))                                       I3                I3                0
     0                                                   0                                                   0                 0                 1];

f = [f1 f2 f3 f4 0]';

G = [-x(2) x(1)  x(4)  -x(3);
     -x(3) -x(4) x(1)  x(2);
     -x(4) x(3)  -x(2) x(1)];

% E = [-x(2) x(1)  -x(4) x(3);
%      -x(3) x(4)  x(1)  -x(2);
%      -x(4) -x(3) x(2)  x(1)];
% 
% A = E*G'

gg = [G'*[x(5);x(6);x(7)]/2 ; M\f];
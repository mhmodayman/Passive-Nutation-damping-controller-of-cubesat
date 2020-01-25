function y = angls(x, ann) 

G = [-x(2) x(1)  x(4)  -x(3);
     -x(3) -x(4) x(1)  x(2);
     -x(4) x(3)  -x(2) x(1)];

E = [-x(2) x(1)  -x(4) x(3);
     -x(3) x(4)  x(1)  -x(2);
     -x(4) -x(3) x(2)  x(1)];

A = E*G';

y = A*ann;
function[H0] = drncos2(t);

% Right ascension of ascending node
omega01 = 272.847*pi/180;     % SV #4 (see Appendix B)
omega02 = 332.847*pi/180;     % SV #7
omega03 = 32.847*pi/180;      % SV #11
omega04 = 212.847*pi/180;     % SV #23

% Mean anomaly
theta01 = 41.806*pi/180;
theta02 = 309.976*pi/180;
theta03 = 339.666*pi/180;
theta04 = 105.206*pi/180;     % END OF PARAMETERS DECLARATION

OmegaRate = 2*pi/86164;
ThetaRate = 2*pi/43082;       % RATE EQUATIONS

omega1 = omega01 - OmegaRate*t;
omega2 = omega02 - OmegaRate*t;
omega3 = omega03 - OmegaRate*t;
omega4 = omega04 - OmegaRate*t;

theta1 = theta01 + ThetaRate*t;
theta2 = theta02 + ThetaRate*t;
theta3 = theta03 + ThetaRate*t;
theta4 = theta04 + ThetaRate*t;   % COMPUTING OMEGA AND THETA

R = 26559800.0;          % ORBIT RADIUS
ERadius = 6380000.0;     % EARTH RADIUS
CosIncl = 0.573576;      % COS 55 DEGREES
SinIncl = 0.819152;      % SIN 55 DEGREES

%  FIRST SATELLITE (SV #4)
x1 = R*(cos(theta1)*sin(omega1)+sin(theta1)*cos(omega1)*CosIncl);
y1 = R*sin(theta1)*SinIncl;
z1 = R*(cos(theta1)*cos(omega1) - sin(theta1)*sin(omega1)*CosIncl)-ERadius;
xyz1 = sqrt (x1^2 + y1^2 + z1^2);


%  SECOND SATELLITE (SV #7)
x2 = R*(cos(theta2)*sin(omega2)+sin(theta2)*cos(omega2)*CosIncl);
y2 = R*sin(theta2)*SinIncl;
z2 = R*(cos(theta2)*cos(omega2) - sin(theta2)*sin(omega2)*CosIncl)-ERadius;
xyz2 = sqrt (x2^2 + y2^2 + z2^2);


%  THIRD SATELLITE (SV #11)
x3 = R*(cos(theta3)*sin(omega3)+sin(theta3)*cos(omega3)*CosIncl);
y3 = R*sin(theta3)*SinIncl;
z3 = R*(cos(theta3)*cos(omega3) - sin(theta3)*sin(omega3)*CosIncl)-ERadius;
xyz3 = sqrt (x3^2 + y3^2 + z3^2);


%  FOURTH SATELLITE (SV #23)
x4 = R*(cos(theta4)*sin(omega4)+sin(theta4)*cos(omega4)*CosIncl);
y4 = R*sin(theta4)*SinIncl;
z4 = R*(cos(theta4)*cos(omega4) - sin(theta4)*sin(omega4)*CosIncl)-ERadius;
xyz4 = sqrt (x4^2 + y4^2 + z4^2);


H0=[-x1/xyz1 -y1/xyz1 -z1/xyz1;
    -x2/xyz2 -y2/xyz2 -z2/xyz2;
    -x3/xyz3 -y3/xyz3 -z3/xyz3;
    -x4/xyz4 -y4/xyz4 -z4/xyz4];

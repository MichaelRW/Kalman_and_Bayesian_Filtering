[4:5] Equation
: Right ascension of ascending node
omega01 = 332.847*pi/180     : SV #7 (see Appendix B)
omega02 = 32.847*pi/180      : SV #10
omega03 = 152.847*pi/180     : SV #17
omega04 = 212.847*pi/180     : SV #24

: Mean anomaly
theta01 = 309.976*pi/180
theta02 = 11.796*pi/180
theta03 = 197.046*pi/180
theta04 = 135.346*pi/180     : END OF PARAMETERS DECLARATION

OmegaRate = 2*pi/86164
ThetaRate = 2*pi/43082   : RATE EQUATIONS

omega1 = omega01 - OmegaRate*t
omega2 = omega02 - OmegaRate*t
omega3 = omega03 - OmegaRate*t
omega4 = omega04 - OmegaRate*t

theta1 = theta01 + ThetaRate*t
theta2 = theta02 + ThetaRate*t
theta3 = theta03 + ThetaRate*t
theta4 = theta04 + ThetaRate*t   : COMPUTING OMEGA AND THETA

R = 26559800.0          : ORBIT RADIUS
ERadius = 6380000.0     : EARTH RADIUS
CosIncl = 0.573576      : COS 55 DEGREES
SinIncl = 0.819152      : SIN 55 DEGREES

:  FIRST SATELLITE (SV #7)
x1 = R*(cos(theta1)*sin(omega1)+sin(theta1)*cos(omega1)*CosIncl)
y1 = R*sin(theta1)*SinIncl
z1 = R*(cos(theta1)*cos(omega1) - sin(theta1)*sin(omega1)*CosIncl)-ERadius
xyz1 = sqrt (sqr (x1) + sqr (y1) + sqr (z1))

$1$1 = -x1/xyz1
$1$2 = -y1/xyz1
$1$3 = -z1/xyz1
$1$4 = 1
$1$5 = 0

:  SECOND SATELLITE (SV #10)
x2 = R*(cos(theta2)*sin(omega2)+sin(theta2)*cos(omega2)*CosIncl)
y2 = R*sin(theta2)*SinIncl
z2 = R*(cos(theta2)*cos(omega2) - sin(theta2)*sin(omega2)*CosIncl)-ERadius
xyz2 = sqrt (sqr (x2) + sqr (y2) + sqr (z2))

$2$1 = -x2/xyz2
$2$2 = -y2/xyz2
$2$3 = -z2/xyz2
$2$4 = 1
$2$5 = 0

:  THIRD SATELLITE (SV #17)
x3 = R*(cos(theta3)*sin(omega3)+sin(theta3)*cos(omega3)*CosIncl)
y3 = R*sin(theta3)*SinIncl
z3 = R*(cos(theta3)*cos(omega3) - sin(theta3)*sin(omega3)*CosIncl)-ERadius
xyz3 = sqrt (sqr (x3) + sqr (y3) + sqr (z3))

$3$1 = -x3/xyz3
$3$2 = -y3/xyz3
$3$3 = -z3/xyz3
$3$4 = 1
$3$5 = 0

:  FOURTH SATELLITE (SV #24)
x4 = R*(cos(theta4)*sin(omega4)+sin(theta4)*cos(omega4)*CosIncl)
y4 = R*sin(theta4)*SinIncl
z4 = R*(cos(theta4)*cos(omega4) - sin(theta4)*sin(omega4)*CosIncl)-ERadius
xyz4 = sqrt (sqr (x4) + sqr (y4) + sqr (z4))

$4$1 = -x4/xyz4
$4$2 = -y4/xyz4
$4$3 = -z4/xyz4
$4$4 = 1
$4$5 = 0

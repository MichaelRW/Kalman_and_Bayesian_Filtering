clear
%File name is pb05_20.m and this is for Problem 5.20.

%Part (a) is analytic.  System poles (char. roots) are at s=0 and
%s=+ and -j*omega (i.e., on the imaginary axis).

%Part (b):  The F and W matrices need to be worked out analytically.
%Then PHI and Q are computed as follows (see Section 5.3).

omega=.2625161
re=6367253
omegax=omega/sqrt(2)
omegaz=omegax
w=(.001*pi/180)^2
f11=[0 omegaz 0;-omegaz 0 omegax;0 -omegax 0]
F=[f11 eye(3);zeros(3) zeros(3)];
GWGT=[zeros(3) zeros(3);zeros(3) w*eye(3)];
A=[-.5*F .5*GWGT;zeros(6) .5*F'];
B=expm(A);
PHIT=B(7:12,7:12);
PHI=PHIT'
Q=PHI*B(1:6,7:12)
comment='Check PHI and Q for reasonableness. Press ENTER to Continue'
pause

%Part (c): We will use filcovar.m here, so we will declare key
%variables to be global.

global n m s PHI Q H R PZEROM COVAR

n=6
m=3
s=11
H=[1 0 0 0 0 0;0 1 0 0 0 0;1 0 1 0 0 0];
R=[(100/re)^2 0 0;0 (100/re)^2 0;0 0 (pi/(60*180))^2];
PZEROM=[(1000/re)^2 0 0 0 0 0
        0 (1000/re)^2 0 0 0 0
        0 0 (.1*pi/180)^2 0 0 0
        0 0 0 (.02*pi/180)^2 0 0
        0 0 0 0 (.02*pi/180)^2 0
        0 0 0 0 0 (.02*pi/180)^2];

filcovar

var11=zeros(1,11);
var22=zeros(1,11);
var33=zeros(1,11);
var44=zeros(1,11);
var55=zeros(1,11);
var66=zeros(1,11);
for i=1:11
   var11(i)=COVAR(1,6*i-5);
   var22(i)=COVAR(2,6*i-4);
   var33(i)=COVAR(3,6*i-3);
   var44(i)=COVAR(4,6*i-2);
   var55(i)=COVAR(5,6*i-1);
   var66(i)=COVAR(6,6*i);
end
poserr1=re*sqrt(var11);
poserr2=re*sqrt(var22);
psizerr=sqrt(var33);
rmsbiasx=sqrt(var44);
rmsbiasy=sqrt(var55);
rmsbiasz=sqrt(var66);

%Plot the two horizontal rms position errors together.
t=0:.5:5;
plot(t,poserr1,'r',t,poserr2,'g')
title('Press ENTER to Continue')
pause

%Plot the psi z rms error by itself.
plot(t,psizerr)
title('Press ENTER to Continue')
pause

%Plot all three bias rms errors together.
plot(t,rmsbiasx,'r',t,rmsbiasy,'g',t,rmsbiasz,'b')
title('Press ENTER to end Problem 5.20')

%The system appears to be completely observable.

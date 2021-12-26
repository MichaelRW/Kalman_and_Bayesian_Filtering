clear
%This file is for Problem 8.6, Parts (a) and (b)

%Part (a) is analytic and the process and measurement models
%are contained in the code below.
 
%Part (b):  First, enter key parameters.
 
n=2
m=1
s=2001
d=10000
dt=1
PHI=[1 0;0 1];
Q=[0 0;0 0];
R=225;
PZEROM=[100^2 0;0 (.02)^2];
 
%Nominal path is along y-axis, and origin is shifted to -10000 point
%for convenience.  The nominal y trajectory is called ynom and it
%increments in 10m steps beginning at k=0.  H is time variable and
%linearization takes place about ynom.
 
%Now prepare for the main i loop.  We will save the p11, p22, and
%p12 terms of the P matrix separately as 1 x 2001 row vectors.  (If
%these row vectors exceed the size limitation in your version of
%MATLAB, try changing the step size to 20m and just do 1001 steps.)
 
PMINUS=PZEROM;
I=eye(2);
H=zeros(1,2);
gain=zeros(2,1);
PPLUS=zeros(2,2);
p11=zeros(1,2001);
p22=zeros(1,2001);
p12=zeros(1,2001);
posvar=zeros(1,2001);
 
for i=1:2001
   i
   ynom=10*(i-1);
   denom=sqrt(d*d+(d-ynom)*(d-ynom));
   H(1,1)=-(d-ynom)/denom;
   H(1,2)=ynom*H(1,1);
   gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-gain*H)*PMINUS;
   PPLUS=.5*(PPLUS+PPLUS');
   p11(i)=PPLUS(1,1);
   p22(i)=PPLUS(2,2);
   p12(i)=PPLUS(1,2);
   posvar(i)=p11(i)+2*ynom*p12(i)+ynom*ynom*p22(i);
   PMINUS=PHI*PPLUS*PHI'+Q;
end
 
%Now do the plots.  Do bias, scale factor, and position separately.
 
k=0:1:2000;
plot(k,sqrt(p11))
title('Press ENTER to Continue')
pause
 
plot(k,sqrt(p22))
title('Press ENTER to Continue')
pause
 
plot(k,sqrt(posvar))
title('Press ENTER to End Part (b)')

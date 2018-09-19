clear
%File name is pb08_9.m and this is for Problem 8.9.

%We will use the smoothing covariance analysis program smcovar.m.
%First, declare global variables and work out key parameters.

global n m s PHI Q H R PZEROM FILTERP SMOOTHP

s=101
n=2
m=1
dt=20
omega0=.012
c=.0043987
F=[0 1;-omega0^2 -omega0*sqrt(2)]
GWGT=[0 0;0 c]
AA=[-dt*F dt*GWGT;zeros(2) dt*F']
B=expm(AA);
PHIT=B(3:4,3:4);
PHI=PHIT'
Q=PHI*B(1:2,3:4)
H=[1 0]
R=[400]
PZEROM=[900 0;0 (30*omega0)^2]

%Now use smcovar.m

smcovar

%Finally, plot the rms estimation errors for states 1 and 2.
%We will make separate plots for rms11 and rms22.

p11=zeros(1,101);
p22=zeros(1,101);
for i=1:101
   p11(i)=SMOOTHP(1,2*i-1);
   p22(i)=SMOOTHP(2,2*i);
end
rms11=sqrt(p11);
rms22=sqrt(p22);

t=0:20:2000;

plot(t,rms11)
title('Press ENTER to Continue')
pause

plot(t,rms22)
title('Press ENTER to end Problem 8.9')

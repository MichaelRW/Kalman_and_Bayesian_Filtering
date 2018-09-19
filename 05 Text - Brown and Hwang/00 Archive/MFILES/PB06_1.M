clear
%File name is pb06_1.m and this is for Problem 6.1.

%We will run this program using filcovar.m in two stages.  There
%is a switch in parameters in going from the first to second stage.

%The purpose of the first stage is to establish the appropriate
%initial Pzero minus for the second stage.
%First, declare global variables.

global n m s PHI Q H R PZEROM COVAR

%The filter parameters for the first stage are:

c=sqrt(.48)
s=41
n=2
m=1
dt=1
F=[0 1;-.36 -.06];
GWGT=[0 0;0 1];
A=[-dt*F dt*GWGT;zeros(2) dt*F'];
B=expm(A);
PHIT=B(3:4,3:4);
PHI=PHIT';
Q=PHI*B(1:2,3:4);
R=[1];
H=[0 c];
PZEROM=zeros(2);

filcovar

%Now check for steady state condition by plotting p11 and p22.

p11=zeros(1,s);
p22=zeros(1,s);
for i=1:s
   p11(i)=COVAR(1,2*i-1);
   p22(i)=COVAR(2,2*i);
end
t=0:1:(s-1);
plot(t,p11,t,p22)
title('Press ENTER to Continue')
pause

%We will now use the last P matrix in COVAR as PZEROM in
%the prediction run.  Also change the R and s variables.

PZEROM=COVAR(1:2,(2*s-1):(2*s));
R=[1.0e20]
s=21

filcovar

%Now plot the rms prediction error.

ysqerr=zeros(1,s);
for i=1:s
   ysqerr(i)=c*c*COVAR(2,2*i);
end
tt=0:1:(s-1)
plot(tt,sqrt(ysqerr))
title('Press ENTER to end Problem 6.1')

%Note that the prediction rms error at 10 sec is about 1.4 m,
%which is significantly less than the 2 meter rms error that
%would exist without optimal prediction.

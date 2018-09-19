clear
%File name is ex06_1.m and this is for Example 6.1.

%This is a simple prediction covariance analysis exercise,so we do
%not need to generate sample realizations of the process.

%First we need to calculate the model parameters.  We will make R
%very large to approximate worthless measurements.  We will use
%filcovar.m, so declare global variables next.

global n m s PHI Q H R PZEROM COVAR

omega0=.012
c=.0043987
dt=10
n=2
m=1
s=31
F=[0 1;-omega0^2 -(sqrt(2))*omega0]
GWGT=[0 0;0 c]
A=[-dt*F dt*GWGT;zeros(2) dt*F'];
B=expm(A);
PHIT=B(3:4,3:4);
PHI=PHIT'
Q=PHI*B(1:2,3:4)
H=[1 0];
R=[1.0e20];
PZEROM=[0 0;0 0];

%Now call filcovar.  (Note that with large R, this is prediction.)

filcovar

var11=zeros(1,31);
for i=1:31
   var11(i)=COVAR(1,2*i-1);
end
rmsrange=sqrt(var11);
t=0:10:300
plot(t,rmsrange)
title('Press ENTER to end Example 6.1')

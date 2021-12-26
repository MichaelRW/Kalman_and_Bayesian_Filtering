%This file is for Problem 3.10
%Specify model parameters
dt=0.2;
F=[0 1;-100 -10];
G=[0;sqrt(10)];
GWGT=G*G';   %Note that W=1
%Now compute PHI and Q
A=dt*[-F GWGT;zeros(2,2) F'];
B=expm(A);
PHI=B(3:4,3:4)';
Q=PHI*B(1:2,3:4);

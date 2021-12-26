% This file is for Problem 4.6
F=[-2 0;0 -4];
G=[-.5 1.5];
GWGT=[.25 -.75;-.75 2.25];
W=[1];
dt=1;
A=dt*[-F GWGT;zeros(2,2) F'];
B=expm(A);
PHIT=B(3:4,3:4);
PHI=PHIT';
Q=PHI*B(1:2,3:4);
%Now do a one step KF
H=[1 1];
R=[2];
PMINUS=[.0625 -.125;-.125 9/32];
gain=PMINUS*H'*inv(H*PMINUS*H'+R);
IMKH=eye(2,2)-gain*H;
PPLUS=IMKH*PMINUS;
PY=PPLUS(1,1)+2*PPLUS(1,2)+PPLUS(2,2);

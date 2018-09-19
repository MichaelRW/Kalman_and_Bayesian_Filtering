%File name is insphiq.m.  This is a special script file to compute
%the INS PHI and Q for for the 9-state model given by Eq.(10.2.10).

%We will use the method given in Chapter 5, Eqs. (5.3.23 - 5.3.26).

global PHI Q dt omegax omegay wacc wgyro

g=9.8;
re=6.37e6;

%Specify the 9 x 9 F matrix

F=zeros(9,9);
F(1,2)=1;
F(2,3)=-g;
F(3,2)=1/re;
F(3,9)=omegax;
F(4,5)=1;
F(5,6)=-g;
F(6,5)=1/re;
F(6,9)=omegay;
F(7,8)=1;

%Specify the 9 x 9 input power spectral density matrix

GWGT=zeros(9,9);
GWGT(2,2)=wacc;
GWGT(3,3)=wgyro;
GWGT(5,5)=wacc;
GWGT(6,6)=wgyro;
GWGT(8,8)=wacc;
GWGT(9,9)=wgyro;

%Now compute PHI and Q

A=[-dt*F dt*GWGT;zeros(9,9) dt*F'];
B=expm(A);
PHIT=B(10:18,10:18);
PHI=PHIT';
Q=PHI*B(1:9,10:18);

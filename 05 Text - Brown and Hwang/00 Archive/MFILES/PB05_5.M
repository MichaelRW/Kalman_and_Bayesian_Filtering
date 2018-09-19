clear
%File name is pb05_5.m and this is for Problem 5.5.

%The H, R, initial state estimate, and initial Pminus should be
%worked out analytically.  The results are:
%          H=[1 1]
%          R=[2]
%          Initial estimate=[0 0]'
%          Initial Pminus=PZEROM=[1/32 0;0 1/8]
%The PHI and Q matrices are evaluated numerically using the methods
%given in Section 5.3.  Partitioned parts of A are:

F=[0 1;-4 -4];
GWGT=[0 0;0 1];

%dt is unity.  Thus the full A matrix is:

A=[-F GWGT;zeros(2,2) F'];

%The desired PHI and Q are then obtained from:

B=expm(A);
PHIT=B(3:4,3:4);
PHI=PHIT'
Q=PHI*B(1:2,3:4)

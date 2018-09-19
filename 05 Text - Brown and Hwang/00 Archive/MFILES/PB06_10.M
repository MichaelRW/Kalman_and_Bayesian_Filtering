clear
%File name is pb06_10.m and this is for Problem 6.10.

%The F, GWGT, and PZERO matrices can be worked out easily with paper-
%and-pencil methods, so we begin with these as known quantities.

%The Markov component has a time constant of 500 sec (50 steps), so
%we should allow about 200 steps to reach steady-state condition.

%We will be using a modified version of filcovar.m where we only
%save the diagonal elements of P.  The modified version is fildiag.m
%and the end result of fildiag.m is a stacked matrix called PDIAG.

%First, specify the global variables.

global n m s PHI Q H R PZEROM PDIAG

%Next, specify the key parameters and work out PHI and Q.

n=3
m=1
s=200
dt=10
H=[1 0 1];
R=[100];
PZEROM=[400 0 0;0 .16 0;0 0 100];
F=[0 1 0;-.0004 -.02 0;0 0 -.002];
GWGT=[0 0 0;0 .0064 0;0 0 .4];
A=[-dt*F dt*GWGT;zeros(3) dt*F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT'
Q=PHI*B(1:3,4:6)

%Part (a): Run fildiag.m for 200 steps with PZEROM as given above.

fildiag

PARTAP=PDIAG;

%Part (b):  Now make another run with PZEROM set to zero.

PZEROM=zeros(3);

fildiag

PARTBP=PDIAG;

%Save H, PHI, and gain for Problem 6.11.

save prob610 H PHI gain

%Now plot the results and check the steady-state condition.
%We are mainly interested in the estimation error for state 1,
%so we will only plot its error variance.

t=0:10:1990;
plot(t,PARTAP(1,:),'r',t,PARTBP(1,:),'g')
title('Press ENTER to end Problem 6.10')

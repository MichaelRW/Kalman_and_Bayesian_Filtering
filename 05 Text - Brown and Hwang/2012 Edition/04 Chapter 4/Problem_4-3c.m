clear
%This file is for Problem 4.3 part (c).

%Parts (a) and (b) are paper-and-pencilproblems, and the results of
%part (a) are needed for part (c).
%First specify beta, W, and dt.

beta=.2
W=10
dt=.1

%The partitioned parts of the A matrix are (see Sec. 5.3):

F=[0 1 0;0 0 1;0 0 -beta];
GWGT=[0 0 0;0 0 0;0 0 W];

%The full A matrix (including the dt multiplier) is:

A=dt*[-F GWGT;zeros(3,3) F'];

%The desired PHI and Q are then obtained from:

B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT'
Q=PHI*B(1:3,4:6)

%Note that the numerical values of PHI and Q are approximately the
%same as those obtained from the equations given in Problem 5.2.

clear
%File name is pb10_1.m and this is for Problem 10.1.

%This is a covariance analysis problem where we only need to look
%at the p11, p22, and p33 terms.  Therefore, we will use fildiag.m
%for this problem.

%Declare global variables.

global n m s PHI Q H R PZEROM PDIAG

%First, workout key parameters and enter them into the workspace.

n=3
m=1
s=1000
PZEROM=zeros(3)
H=[1 0 1]
R=[100]
dt=1
beta=1/300

%Now calculate PHI and Q

F=[0 1 0;0 0 0;0 0 -beta]
GWGT=[0 0 0;0 .13889 0;0 0 2*10000*beta]
A=[-dt*F dt*GWGT;zeros(3) dt*F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT'
Q=PHI*B(1:3,4:6)

%Now call fildiag.m.  (Note that fildiag.m will produce a 3 x 1000
%matrix called PDIAG.  If the matrix size exceeds the limits of your
%version of MATLAB, we suggest that you modify fildiag.m and define
%the p11, p22, and p33 rows of PDIAG as separate row vectors with
%different variable names.  Then modify the subsequent code in 
%the main program accordingly.)

fildiag

%Plot p11 and p33 together for comparison.

t=0:1:999;
plot(t,PDIAG(1,:),t,PDIAG(3,:))
title('Press ENTER to Continue')
pause

%Plot p22 separately.  (Units are different from p11 and p33.)

plot(t,PDIAG(2,:))
title('Press ENTER to end Problem 10.1')

clear
%File name is pb07_9.m and this is for Problem 7.9.

%Part (a):  This is a paper-and-pencil problem.  The solution is:
%           p11=p22=100*sqrt(2), p12=100; and 
%           eigenvalues=-(1/sqrt(2)) + and - j*(1/sqrt(2))

%Part (b):  Enter key parameters into workspace in anticipation of
%using filcovar.m.  We will first try 25 steps and then run the filter
%again for 50 steps and compare to see if steady-state is achieved.

%Declare global variables.

global n m s PHI Q H R PZEROM COVAR

PHI=[1 .5;0 1]
Q=[100/24 100/8;100/8 50]
R=[200]
H=[1 0]
PZEROM=zeros(2)
s=25
n=2
m=1

%Now run filcovar.m and make a note of final P matrix.

filcovar
COVAR(1:2,49:50)
comment='Note p11 and p22.  Press ENTER to continue.'
pause

%Rerun filcovar.m for 50 steps after re-intializing s.

s=50
filcovar
COVAR(1:2,99:100)
comment='p11 and p22 are same as before. Press ENTER to continue.'
pause

%Filter appears to have reached steady-state after 50 steps.
%Using the gain at step 50, compute the eigenvalues.

A=PHI-gain*H*PHI;
b=eig(A)

%Conclusion:  The characteristic poles of the discrete system are
%similar to those of the continuous system when one considers the
%relatioship between the s-plane and the z-plane.

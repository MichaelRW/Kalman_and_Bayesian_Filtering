clear
%File name is pb05_7.m and this file is for Problem 5.7.
                                                                
%We will use the covariance analysis program filcovar.m here.
%Declare global variables common to filcovar.m and main program.

global n m s PHI Q H R PZEROM COVAR

%Part (a):  Your general program (or filcovar.m) should yield the
%same results as obtained in problem 5.6.

%Part (b):  Some of the parameters in the model are to be worked 
%out analytically (see Problem 5.2).  The full set of parameters
%to be entered into the workspace before using filcovar.m are:

n=2
m=1
s=51
PHI=[1 1;0 1]
Q=[10/3 10/2;10/2 10]
H=[1 0]
R=[225]
PZEROM=[0 0;0 0]

%Now run filcovar.m

filcovar

%Plot the p11 and p22 sequences and look for a steady-state condition.

p11=zeros(1,51);
p22=zeros(1,51);
for i=1:51
   p11(i)=COVAR(1,(2*i-1));
   p22(i)=COVAR(2,2*i);
end
t=0:1:50;
plot(t,p11,'r',t,p22,'g')
title('Press ENTER to end Problem 5.7')

%It appears that the system is observable.  Otherwise, one or
%the other of the error covariances would increase without
%bound, because both x1 and x2 are nonstationary.

clear
%File name is pb06_12.m and this is for Problem 6.12.

%This is a simple error covariance analysis example.  We do not need the
%off-diagonal terms, so we will use fildiag.m (rather than filcovar.m)
%and run the filter for 100 steps and observe the trends of p11 and p22.

%First, declare global variables.

global n m s PHI Q H R PZEROM PDIAG

%n is the dimension of the state vector.
%m is the dimension of the measurement vector.
%s is the number of steps (including k=0).

%Next, enter the key parameters.

n=2
m=1
s=100
PHI=[1 1;0 1]
Q=.1*[1/3 1/2;1/2 1]
H=[0 1]
R=[.01]
PZEROM=zeros(2)

%Run fildiag.m and plot the results.

fildiag

t=0:1:(s-1);
plot(t,PDIAG(1,:),'r',t,PDIAG(2,:),'g')
title('Press ENTER to end Problem 6.12')

%Clearly, p11 is ramping off linearly, so there is apparent
%"divergence".  However, the filter is doing the best that it
%can in a poor observability situation.

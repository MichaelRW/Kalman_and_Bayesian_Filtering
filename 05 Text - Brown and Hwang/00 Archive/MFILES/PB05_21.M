clear
%File name is pb05_21.m and this is for Problem 5.21.

%This is a covariance analysis exercise, so we will use filcovar.m
%for the 1st part.  First, look at measurement model 1.  We need
%to declare the global variables.

global n m s PHI Q H R PZEROM COVAR

n=1
m=1
s=201
PHI=1
Q=4
H=.5
R=1
PZEROM=100

filcovar

result1=COVAR;
save prob521 result1
k=0:1:200;
plot(k,result1)
title('Press ENTER to Continue')
pause

%The H matrix in the second measurement model varies with k.  File
%filcovar.m is set up for constant H, so it cannot be used directly.
%The easiest way to accommodate variable H here is to create a new
%file, say covar521.m, where H varies with i within the "for" loop.
%File covar521.m is provided for this purpose.
%Clear and declare a new set of global variables for 2nd part.

clear
global n m s PHI Q R PZEROM COVAR
n=1
m=1
s=201
PHI=1
Q=4
R=1
PZEROM=100

covar521

result2=COVAR;
k=0:1:200;
load prob521
plot(k,result1,'r',k,result2,'g')
title('Press ENTER to end Problem 5.21')
      
%The coupling between the observable and the state becomes poor
%around k=70, and thus the large error variance at this point.

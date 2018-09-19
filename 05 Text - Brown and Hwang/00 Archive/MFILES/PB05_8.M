clear
%File name is pb05_8.m and this is for Problem 5.8 part (b).

%Part (a) is a paper-and-pencil problem.  The states are decoupled,
%so it is easy to work out the filter parameters in general form.

%Part (b):  We will be using the filter covariace analysis program
%filcovar.m, so we will declare the filter parameters and the
%filcovar output matrix as global variables.

global n m s PHI Q H R PZEROM COVAR

%n is the state vector dimension.
%m is the dimension of the mesurement vector.
%s is the number of steps (including the step at k=0).

PHI=[exp(-.1*1) 0;0 exp(-1*1)]
Q=[9*(1-exp(-2*.1*1)) 0;0 1*(1-exp(-2*1*1))]
H=[1 1]
R=[0]
PZEROM=[9 0;0 1]

n=2
m=1
s=51

%Now run filcovar.m

filcovar

%Examine the p11 and p22 sequences for reasonableness.  These
%can be easily extracted from the stacked matrix COVAR as follows:

p11=zeros(1,51);
p22=zeros(1,51);
for i=1:51
   p11(i)=COVAR(1,(2*i-1));
   p22(i)=COVAR(2,2*i);
end

%p11 and p22 can now be plotted easily.

t=0:1:50;
plot(t,p11,'r',t,p22,'g')
title('Press ENTER to end Problem 5.8')

%Note that the p11 and p22 sequences are identical.  Is this 
%reasonable?  Explain.

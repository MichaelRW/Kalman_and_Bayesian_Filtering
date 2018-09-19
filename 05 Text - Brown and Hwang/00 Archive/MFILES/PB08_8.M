clear
%File name is pb08_8.m and this is for Problem 8.8.

%The general smoothing error covariance analysis program smcovar.m
%will be used.  Declare global variables and specify parameters.

global n m s PHI Q H R PZEROM FILTERP SMOOTHP

n=1
m=1
s=51
PHI=exp(-.02)
Q=1-exp(-.04)
H=1
R=1
PZEROM=1

%Now run smcovar.m

smcovar

%Next, plot 51-element vector SMOOTHP.

t=0:.02:1;
plot(t,SMOOTHP,'ow')
title('Press ENTER to end Problem 8.8')

%The estimation accuracy is governed by the correlation structure,
%and it is the same for both positive and negative "lags".  Thus,
%we would expect to have symmetry about the midpoint of the data span.

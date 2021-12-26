clear
%This file is for Problem 3.16.

%We need to do the integration for 100 values of the upper limit t,
%and we will do each integration in a "for" loop using quad8.

%The squared impulsive response of the filter is called gsq329(x),
%and it is in the file gsq329.m.  (The input PSD scale factor is
%included in the function gsq329(x).)

%Now do the integration for 100 points beginning with t=.001.

avexsq=zeros(1,101);
dt=.001
A=10
w0=20*pi
for i=1:100
   i
   a=quad8('gsq329',0,i*dt);
   avexsq(i+1)=a;
end

%Now generate data for the exact response (from Prob. 3.12b).

exact=zeros(1,101);
for i=1:100
   exact(i+1)=(A/(w0^3))*(.5*w0*i*dt-.25*sin(2*w0*i*dt));
end

%Plot the result.

t=0:.001:.1;
plot(t,avexsq,'r',t,exact,'w')
title('Press ENTER to end Problem 3.29')

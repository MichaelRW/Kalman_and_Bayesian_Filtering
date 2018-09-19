clear
%File name is pb09_6a.m and this is for the optimal part of Problem 9.6.

%We separate this problem into three parts: (a) optimal filter, (b) Schmidt-
%Kalman filter, and (c) R-Bumped-Up filter.  Each is in a separate file.

%Working out the error covariance for the optimal filter is routine.
%We will use filcovar.m and save the variance p11 in a mat file.

%First, declare global variables and specify key parameters.

global n m s PHI Q H R PZEROM COVAR

n=2
m=1
s=51
PHI=[1 0;0 exp(-.1)]
Q=[1 0;0 (1-exp(-.2))]
H=[1 1]
R=[1]
PZEROM=[1.618034 0;0 1]

%Now run filcovar.m

filcovar

%Sift out the p11 sequence.

p11opt=zeros(1,51);
for i=1:51
   p11opt(i)=COVAR(1,2*i-1);
end

%The desired error variance will be saved in prob96a.mat

save prob96a p11opt
comment='Go on to the Schmidt-Kalman part.  Call pb09_6b'

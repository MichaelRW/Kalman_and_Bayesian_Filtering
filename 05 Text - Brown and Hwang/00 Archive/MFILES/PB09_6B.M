clear
%File name is pb09_6b.m and this is for the Schmidt-Kalman part
%of Problem 9.6.

%We will use the general sch_kal.m analysis program here.  The first
%step is to declare global variables and specify key parameters.

global n r m s PHIX PHIY QX QY H J R PX0M PXY0M PY0M COVARX

n=1
r=1
m=1
s=51
PHIX=1
PHIY=exp(-.1)
QX=1
QY=1-exp(-.2)
H=1
J=1
R=1
PX0M=1.618034
PXY0M=0
PY0M=1

%Now use sch_kal.m

sch_kal

%Sift out the desired error variance sequence (just COVARX in this
%problem), and save in file prob96b.mat

p11sk=COVARX;
save prob96b p11sk
comment='Go on to R-Bumped-Up part.  Call pb09_6c'

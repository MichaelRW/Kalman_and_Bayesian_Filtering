%File name is ex06_5.m and this is for Example 6.5.

%This is a continuation of Example 6.1, so the simplest way to get the
%comparison plot is to reuse ex06_1.m and add a few more statements.

ex06_1

%We need to calculate PHISUB, PHISA, and QSA for 21 steps.  We will
%do this in a "for" loop and stack the results side-by-side in PSUB.

COVARX=[900 0;0 .1296]
PSUB=zeros(2,42);
sqerrsub=zeros(1,21);
rmssub=zeros(1,21);
for i=1:20
   i
   ASA=[-dt*i*F dt*i*GWGT;zeros(2) dt*i*F'];
   BSA=expm(ASA);
   PHISAT=BSA(3:4,3:4);
   PHISA=PHISAT';
   QSA=PHISA*BSA(1:2,3:4);
   PHISUB=[1 dt*i;0 1];
   PHIDIFF=PHISA-PHISUB;
   PSUB(1:2,(2*i+1):(2*i+2))=PHIDIFF*COVARX*PHIDIFF'+QSA;
   sqerrsub(i+1)=PSUB(1,2*i+1);
   rmssub=sqrt(sqerrsub);
end

%Now plot the optimal and suboptimal rms errors together.

tt=0:10:200;
plot(tt,rmsrange(1:21),'r',tt,rmssub,'g')
title('Press ENTER to end Example 6.5')

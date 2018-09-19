clear
%File name is pb09_5.m and this is for Problem 9.5.

%Part (a):  This part is analytic.  Answer is given in problem statement.

%Part (b):  This is covariance analysis using delayed-state algorithm.
%We only need to save p11, p22, and p33, so this program will be
%similar to fildiag.m.  Final result is stacked matrix PDIAG.

%First work out and put key parameters into the workspace.  The
%notation is the same as in fildiag.m.

n=3
m=1
s=101
g=32.2
re=2.09e7
omega0=sqrt(g/re)
dt=120
PZEROM=zeros(3);

%Compute PHI and Q

F=[0 1 0;0 0 -g;0 1/re 0]
GWGT=[0 0 0;0 .003 0;0 0 2.35e-11]
A=[-dt*F dt*GWGT;zeros(n) dt*F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT';
Q=PHI*B(1:3,4:6);

%We will use two different values of R.  H and J are constant.

R1=1.0e20
R2=90000
H=[1 0 0]
J=[-1 0 0]

%Now set up the "for" loop for s steps.

PDIAG=zeros(n,s);
PMINUS=PZEROM;
PPLUS=PZEROM;
for i=1:s
   i
   if i<=43
      R=R1;
   else
      R=R2;
   end
   L=H*PMINUS*H'+R+J*PPLUS*PHI'*H'+H*PHI*PPLUS*J'+J*PPLUS*J';
   gain=(PMINUS*H'+PHI*PPLUS*J')*inv(L);
   PPLUS=PMINUS-gain*L*gain';
   PPLUS=.5*(PPLUS+PPLUS');
   for j=1:n
      PDIAG(j,i)=PPLUS(j,j);
   end
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Now plot the results.  Look at p11, p22, and p33 separately.

k=0:1:(s-1);

plot(k,PDIAG(1,:))
title('Press ENTER to Continue')
pause

plot(k,PDIAG(2,:))
title('Press ENTER to Continue')
pause

plot(k,PDIAG(3,:))
title('Press ENTER to Continue')
pause

%Part (c):  Use the PHI, gain, H, and J from Part (b) for the
%eigenvalue calculation.

AA=PHI-gain*H*PHI-gain*J;
chroots=eig(AA)
comment='Roots are as expected. This ends Problem 9.5'

%Part(d):  This exercise is accomplished in a routine manner by
%simply editting R2 to the desired numerical values.

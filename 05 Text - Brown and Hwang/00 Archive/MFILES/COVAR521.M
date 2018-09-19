%Purpose is to compute the covariance with variable H for Pb. 5.21.

%File name is covar521.m.
%Declare global variables.

global n m s PHI Q R PZEROM COVAR

%n is the state dimension.
%m is the dimension of the measurement vector.
%s is the number of steps (including the step at k=0).

%First form the general identity matrix.

I=zeros(n,n);
for i=1:n
   I(i,i)=1;
end

COVAR=zeros(n,n*s);
PMINUS=PZEROM;
for i=1:s
   H=cos(1+i/120);
   gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-gain*H)*PMINUS;
   PPLUS=.5*(PPLUS+PPLUS');
   COVAR(1:n,(n*(i-1)+1):n*i)=PPLUS;
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Note that the first (n x n) block of COVAR corresponds to PPLUS
%at k=0, the second block is for k=1, etc.

%Also note that PPLUS is symmetrized after each update.

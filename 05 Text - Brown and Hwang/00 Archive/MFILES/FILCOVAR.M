%The purpose of this program is to compute the error covariance
%sequence for a discrete Kalman filter with constant parameters.
%The final result is a horizontal stacked matrix called COVAR.

%File name is filcovar.m

%Key parameters n, m, s, PHI, Q, H, R, PZEROM, and the output of
%filcovar.m ,COVAR, are declared as global variables.

global n m s PHI Q H R PZEROM COVAR

%n is the state vector dimension.
%m is the dimension of the measurement vector.
%s is the number of steps (including the step at k=0).

%First form the general identity matrix (no ambiguity).

I=zeros(n,n);
for i=1:n
   I(i,i)=1;
end

%Set up the "for" loop for s steps.

COVAR=zeros(n,n*s);
PMINUS=PZEROM;
for i=1:s
   gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-gain*H)*PMINUS;
   %Symmetrize PPLUS
   PPLUS=.5*(PPLUS+PPLUS');
   COVAR(1:n,(n*(i-1)+1):n*i)=PPLUS;
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Note that the first (n x n) block of COVAR corresponds to PPLUS
%at k=0, the second block is for k=1, etc.

%Also note that PPLUS is symmetrized after each update.

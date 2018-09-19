%Purpose is to compute Kalman filter P and save the diagonal terms.
%The final result is a horizontal stacked matrix PDIAG
%containing the diagonal terms of P as column vectors.

%File name is fildiag.m

%First declare the global variables.

global n m s PHI Q H R PZEROM PDIAG

%n is the state dimension.
%m is the dimension of the measurement vector.
%s is the number of steps (including the step at k=0).

%First form the general identity matrix.

I=zeros(n,n);
for i=1:n
   I(i,i)=1;
end

PDIAG=zeros(n,s);
PMINUS=PZEROM
for i=1:s
   gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-gain*H)*PMINUS;
   PPLUS=.5*(PPLUS+PPLUS');
   for j=1:n
      PDIAG(j,i)=PPLUS(j,j);
   end
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Note that the first (n x 1) column of PDIAG corresponds to PPLUS
%at k=0, the second column is for k=1, etc.

%Also note that PPLUS is symmetrized after each update.

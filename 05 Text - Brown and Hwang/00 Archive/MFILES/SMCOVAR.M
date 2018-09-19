%This is general program for finding both the discrete filter and
%smoothing error covariances for a fixed interval.  The Rauch, Tung,
%Striebel algorithm is used.  The stacked smoothed P is SMOOTHP.
                                                               
%File name is smcovar.m

%Declare global variables.
                                      
global n m s PHI Q H R PZEROM FILTERP SMOOTHP

%s=Number of steps (including k=0).
%n=Dimension of state vector.
%m=Dimension of measurement vector.

%First, do the forward sweep and save both the prior and updated P
%sequences in stacked matrices.

%Form the general identity matrix I.

I=zeros(n,n);
for i=1:n
   I(i,i)=1;
end

PRIORP=zeros(n,n*s);
FILTERP=zeros(n,n*s);
SMOOTHP=zeros(n,n*s);
%PM is the same as PMINUS in previous filcovar and fildiag programs.
PM=PZEROM;
for i=1:s
   i
   filgain=PM*H'*inv(H*PM*H'+R);
   PPLUS=(I-filgain*H)*PM;
   PPLUS=.5*(PPLUS+PPLUS);
   PRIORP(1:n,(n*(i-1)+1):(n*i))=PM;
   FILTERP(1:n,(n*(i-1)+1):(n*i))=PPLUS;
   PM=PHI*PPLUS*PHI'+Q;
end

%Now do the backward sweep.  The end result is the stacked matrix
%called SMOOTHP.  The n x n component matrices in the stack are
%in ascending order in time (just as with the filter P matrices.)

A=zeros(n);
SMOOTHP(1:n,(n*(s-1)+1):(n*s))=FILTERP(1:n,(n*(s-1)+1):(n*s));
for i=1:(s-1)
   j=n*(s-i)-n+1
   k=n*(s-i);
   A=FILTERP(1:n,j:k)*PHI'*inv(PRIORP(1:n,(j+n):(k+n)));
   SMOOTHP(1:n,j:k)=FILTERP(1:n,j:k)+A*(SMOOTHP(1:n,(j+n):(k+n))-...
                    PRIORP(1:n,(j+n):(k+n)))*A';
end

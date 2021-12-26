%This file is for Problem 6.4
%Specify parameters for RTS algorithm
sigma=1.0;
omega0=.1;
bsq=2*(sqrt(2))*(omega0^3);
dt=1.0;
b=sqrt(bsq);
F=[0 1;-omega0^2 -(sqrt(2))*omega0];
G=[0;b*sigma];
W=1.0;
GWGT=G*W*G';
A=dt*[-F GWGT;zeros(2) F'];
B=expm(A);
PHI=B(3:4,3:4)';
Q=PHI*B(1:2,3:4);
H=[1 0];
R=.25;
PMINUS0=[sigma^2 0;0 (omega0*sigma)^2];
N=50;   %number of steps not including step at t=0
%Now do the RTS smoother covariance
%Save the elements of P matrices as 50 point time vectors
pminus11=zeros(N,1);
pminus12=zeros(N,1);
pminus22=zeros(N,1);
pplus11=zeros(N,1);
pplus12=zeros(N,1);
pplus22=zeros(N,1);
%Do the t=0 solution separately
gain0=PMINUS0*H'*inv(H*PMINUS0*H'+R);
IMKH=(eye(2)-gain0*H);
PPLUS0=IMKH*PMINUS0;
PPLUS=PPLUS0;
%Now do forward filter pass
for i=1:N
  PMINUS=PHI*PPLUS*PHI'+Q;
  gain=PMINUS*H'*inv(H*PMINUS*H'+R);
  IMKH=(eye(2)-gain*H);
  PPLUS=IMKH*PMINUS;
  %Save the elements of P matrices
  pminus11(i)=PMINUS(1,1);
  pminus12(i)=PMINUS(1,2);
  pminus22(i)=PMINUS(2,2);
  pplus11(i)=PPLUS(1,1);
  pplus12(i)=PPLUS(1,2);
  pplus22(i)=PPLUS(2,2);
end
%Next,do the backward pass to get smoother solution
%PPLUS at k=50 is also the smoother solution at k=50
PSMTH=PPLUS;   %Backward initial condition

psmth11=zeros(N,1);
psmth11(N)=PSMTH(1,1);
for i=1:(N-1)
  PPLUS=[pplus11(N-i) pplus12(N-i);pplus12(N-i) pplus22(N-i)];
  PMINUS=[pminus11(N-i+1) pminus12(N-i+1);pminus12(N-i+1) ...
         pminus22(N-i+1)];
  ASMTH=PPLUS*PHI'*inv(PMINUS);
  PSMTH=PPLUS+ASMTH*(PSMTH-PMINUS)*ASMTH';
  %Save state 1 smoother variance
  psmth11(N-i)=PSMTH(1,1);
end
%Smoother variance at t=0 is the same as that at t=50
filrms=sqrt([PPLUS0(1,1);pplus11]);
smrms=sqrt([pplus11(50);psmth11]);
%plot filrms and smrms



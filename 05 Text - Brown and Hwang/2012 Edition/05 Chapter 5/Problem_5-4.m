%This file is for Problem 5.4
%Filter model is 2-state. Specify parameters
W=.1;   %psd of input white noise
dt=1.0;
PHI=[1 dt;0 1];
q11=(W/3)*(dt^3);
q12=(W/2)*(dt^2);
q21=q12;
q22=W*dt;
Q=[q11 q12;q21 q22];
H=[0 1];
R=.01;
n=100;
PPLUS=zeros(2,2);
PMINUS=zeros(2,2);
p11=zeros(n,1);
p22=zeros(n,1);
%Now do the covariance run
for i=1:n
  PMINUS=PHI*PPLUS*PHI'+Q;
  GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
  IMKH=(eye(2)-GAIN*H);
  PPLUS=IMKH*PMINUS;
  %Save p11 and p22
  p11(i)=PPLUS(1,1);
  p22(i)=PPLUS(2,2);
end

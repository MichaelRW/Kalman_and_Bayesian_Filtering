%This file is for Problem 5.1
%This is (a) part. Specify parameters
PHI=[1 0;0 1];
Q=[100 0;0 0];
PMINUS0=[1.0e6 0;0 1.0e6];
R=16;
h0=[1 1];
%Do t=0 step separately
gain0=PMINUS0*h0'*inv(h0*PMINUS0*h0'+R);
IMKH=eye(2)-gain0*h0;
rmserr1=zeros(600,1);
rmserr2=zeros(600,1);
PPLUS0=IMKH*PMINUS0;
PMINUS=PHI*PPLUS0*PHI'+Q;
for k=1:600
  h=[cos(2*pi*k/600) 1];
  gain=PMINUS*h'*inv(h*PMINUS*h'+R);
  IMKH=eye(2)-gain*h;
  PPLUS=IMKH*PMINUS;
  %Save rmserrors
  rmserr1(k)=sqrt(PPLUS(1,1));
  rmserr2(k)=sqrt(PPLUS(2,2));
  PMINUS=PHI*PPLUS*PHI'+Q;
end
%Define 601-point rms errors for plotting
rms1601=[sqrt(PPLUS0(1,1));rmserr1];
rms2601=[sqrt(PPLUS0(2,2));rmserr2];



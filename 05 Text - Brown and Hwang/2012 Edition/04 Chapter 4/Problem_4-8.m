%This file is for Problem 4.8
%Specify parameters
sigmas=3.0;
sigman=1.0;
betas=.1;
betan=1.0;
dt=1.0;
PHI=[exp(-betas*dt) 0;0 exp(-betan*dt)];
Q=[(sigmas^2)*(1-exp(-2*betas*dt)) 0;...
   0 (sigman^2)*(1-exp(-2*betan*dt))];
H=[1 1];
R=0.0;
PMINUS0=[sigmas^2 0;0 sigman^2];
%Now do an update at t=0
gain0=PMINUS0*H'*inv(H*PMINUS0*H'+R);
PPLUS0=(eye(2)-gain0*H)*PMINUS0;
PMINUS=PHI*PPLUS0*PHI'+Q;    %To be used at next step
%Now do 50 more steps beginning at t=1
p11=zeros(50,1);
p22=zeros(50,1);
for i=1:50
  gain=PMINUS*H'*inv(H*PMINUS*H'+R);
  IMKH=eye(2)-gain*H;
  PPLUS=IMKH*PMINUS*IMKH'+gain*R*gain';
  %Now save the error variances
  p11(i)=PPLUS(1,1);
  p22(i)=PPLUS(2,2);
  PMINUS=PHI*PPLUS*PHI'+Q;
end















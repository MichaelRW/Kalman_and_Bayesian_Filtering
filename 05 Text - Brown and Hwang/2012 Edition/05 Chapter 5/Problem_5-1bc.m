%This file is for Problem 5.1
%This is parts b and c.  Specify parameters
PHI=[1 0;0 1];
Q=[100 0;0 1.0e-8];
QINV=inv(Q);
PMINUS0=[1.0e6 0;0 1.0e6];
INMINUS0=inv(PMINUS0);
R=16;
RINV=1/16;
h0=[1 1];
%Do t=0 step separately
INPPLUS0=INMINUS0+h0'*RINV*h0;
PPLUS0=inv(INPPLUS0);
M0=(inv(PHI'))*INPPLUS0*inv(PHI);
INMINUS=M0-M0*(inv(M0+QINV))*M0';
rmserr1=zeros(600,1);
rmserr2=zeros(600,1);
PHIINV=inv(PHI);
PHITINV=inv(PHI');
for k=1:600
  h=[cos(2*pi*k/600) 1];
  INPPLUS=INMINUS+h'*RINV*h;
  PPLUS=inv(INPPLUS);
  %Save rms errors
  rmserr1(k)=sqrt(PPLUS(1,1));
  rmserr2(k)=sqrt(PPLUS(2,2));
  %Now project ahead and get INMINUS
  M=PHITINV*INPPLUS*PHIINV;
  INMINUS=M-M*(inv(M+QINV))*M';
end
%Finally, define 601-point vectors for
%plotting and comparison with part (a).
rms1601=[sqrt(PPLUS0(1,1));rmserr1];
rmserr2=[sqrt(PPLUS0(2,2));rmserr2];
   
  

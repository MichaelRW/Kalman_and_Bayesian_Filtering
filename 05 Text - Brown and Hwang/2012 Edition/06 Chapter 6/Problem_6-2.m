%This file is for Problem 6.2
%Specify parameters
sigmax=30.0;
betax=1.0;
dt=.02;
h=1.0;
r=1.0;
phi=exp(-betax*dt);
q=(sigmax^2)*(1-exp(-2*betax*dt));
%Do filter solution separately at t=0
pminus0=sigmax^2;
gain0=pminus0*h/(h*pminus0*h+r);
imkh=1-gain0*h;
pplus0=imkh*pminus0;
%Define vectors for forward sweep. Use Sec. 6.2 notation
N=50;
pminus=zeros(51,1);
pplus=zeros(50,1);
pminus(1)=phi*pplus0*phi+q;
for k=1:N
  gain=pminus(k)*h/(h*pminus(k)*h +r);
  imkh=1-gain*h;
  pplus(k)=imkh*pminus(k);
  pminus(k+1)=phi*pplus(k)*phi+q;   %Note final pminus is not used
end
%Now do backward sweep to get psmooth
%pplus(50) is smoother solution at k=50
psmooth=zeros(50,1);
psmooth(50)=pplus(50);
A=zeros(N,1);
for i=1:(N-1)
  A(N-i)=pplus(N-i)*phi/pminus(N+1-i);
  psmooth(N-i)=pplus(N-i)+A(N-i)*(psmooth(N-i+1)-pminus(N-i+1))*A(N-i);
end
%Now complete the smoother solution at t=0
A0=pplus0*phi/pminus(1);
psmooth0=pplus0+A0*(psmooth(1)-pminus(1))*A0;

  

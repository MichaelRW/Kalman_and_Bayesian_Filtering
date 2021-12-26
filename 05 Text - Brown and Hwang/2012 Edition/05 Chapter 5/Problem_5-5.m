%This file is for Problem 5.5
%Extend incorrect model with r=1
%Specify parameters
sigmax=1.0;
betax=.1;
dt=1.0;
incr=1.0;
truer=2.0;
phi=exp(-betax*dt);
q=(sigmax^2)*(1-exp(-2*betax*dt));
h=1.0;
n=100;
%Solution at t=0 is from Example 5.3
incgain=.5;
incp=.5;
pplus=incp;
subgain=zeros(n,1);
%Now extend incorrect solution for 100 steps
for i=1:100
  pminus=phi*pplus*phi+q;
  gain=pminus*h/(h*pminus*h+incr);
  imkh=1-gain*h;
  pplus=imkh*pminus;
  %save gain sequence in subgain
  subgain(i)=gain;
end
%Now get realistic p by recycling sub gains thru truth model
%Solution at t=0 is from paper and pencil methods
preal0=.75;
preal=zeros(n,1);
pplus=preal0;
for j=1:100
  pminus=phi*pplus*phi+q;
  imkh=1-h*subgain(j);
  pplus=imkh*pminus*imkh+subgain(j)*truer*subgain(j);
  preal(j)=pplus;
end
%Finally, get optimal solution using true r=2
pplus=2/3;
popt=zeros(n,1);
for k=1:100
  pminus=phi*pplus*phi+q;
  gain=pminus*h/(h*pminus*h+truer);
  imkh=1-gain*h;
  pplus=imkh*pminus;
  %save optimal p
  popt(k)=pplus;
end
t=1:100
X=[preal popt];
plot(X,t)




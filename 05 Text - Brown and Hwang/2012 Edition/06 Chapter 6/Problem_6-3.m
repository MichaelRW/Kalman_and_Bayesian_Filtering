%This file is for Problem 6.3
%Fraser-Potter smoothing problem --chap 6
%specify parameters

dt=.02;
sigma=1.0;
sigmasq=sigma^2;
beta=1.0;
phi=exp(-beta*dt);
hk=1.0;
rk=1.0;
qk=sigmasq*(1-exp(-2*beta*dt));
N=51;
%First run the part1 forward filter 26 steps
pplus1=zeros(26,1);
pminus=1.0;
for i=1:26
  gain1=pminus*hk/(hk*pminus*hk+rk);
  imkh1=1-gain1*hk;
  pplus1(i)=imkh1*pminus;
  pminus=phi*pplus1(i)*phi+qk;
end
%Now run the part2 backward filter 25 steps
pplus2=zeros(25,1);
for i=1:25
  if i==1
    pplus2(1)=1.0;
    pminus=phi*pplus2(1)*phi+qk;
  else
    gain2=pminus*hk/(hk*pminus*hk+rk);
    imkh2=1-gain2*hk;
    pplus2(i)=imkh2*pminus;
    pminus=phi*pplus2(i)*phi+qk;
  end
end
%Final step is to form smoother covariance at k=25
psmooth=pplus1(26)*pminus/(pplus1(26)+pminus);


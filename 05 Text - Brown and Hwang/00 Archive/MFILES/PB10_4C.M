clear
%File name is pb10_4c.m and this is for Problem 10.4 , Part (c).

%Enter key parameters

n=2
m=1
s=2001
d=10000
dt=1
PHI=[1 0;0 1];
Q=[0 0;0 0];
PZEROM=[100^2 0;0 (.02)^2];

%Nominal path is along y-axis, and origin is shifted to -10000 point
%for convenience.  The nominal y trajectory is called ynom and it
%increments in 10m steps beginning at k=0.  H is time variable and
%linearization takes place about ynom.  R switches in middle of run.

%Now prepare for the main i loop.  We will save the p11, p22, and
%p12 terms of the P matrix separately as 1 x 2001 row vectors.  (If
%these row vectors exceed the size limitation in your version of
%MATLAB, try changing the step size to 20m and just do 1001 steps.)

PMINUS=PZEROM;
I=eye(2);
H=zeros(1,2);
gain=zeros(2,1);
PPLUS=zeros(2,2);
p11=zeros(1,2001);
p22=zeros(1,2001);
p12=zeros(1,2001);
posvar=zeros(1,2001);

for i=1:2001
   i
   if i>0 & i<=1000
      R=225;
      ynom=10*(i-1);
      denom=sqrt(d*d+(d-ynom)*(d-ynom));
      H(1,1)=-(d-ynom)/denom;
      H(1,2)=ynom*H(1,1);
      gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
      PPLUS=(I-gain*H)*PMINUS;
      PPLUS=.5*(PPLUS+PPLUS');
      p11(i)=PPLUS(1,1);
      p22(i)=PPLUS(2,2);
      p12(i)=PPLUS(1,2);
      posvar(i)=p11(i)+2*ynom*p12(i)+ynom*ynom*p22(i);
      PMINUS=PHI*PPLUS*PHI'+Q;
   end

   if i>1000 & i<=1800
      R=1.0e20;
      ynom=10*(i-1);
      denom=sqrt(d*d+(d-ynom)*(d-ynom));
      H(1,1)=-(d-ynom)/denom;
      H(1,2)=ynom*H(1,1);
      gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
      PPLUS=(I-gain*H)*PMINUS;
      PPLUS=.5*(PPLUS+PPLUS');
      p11(i)=PPLUS(1,1);
      p22(i)=PPLUS(2,2);
      p12(i)=PPLUS(1,2);
      posvar(i)=p11(i)+2*ynom*p12(i)+ynom*ynom*p22(i);
      PMINUS=PHI*PPLUS*PHI'+Q;
   end

   if i>1800 & i<=2001
      R=225;
      ynom=10*(i-1);
      denom=sqrt(d*d+(d-ynom)*(d-ynom));
      H(1,1)=-(d-ynom)/denom;
      H(1,2)=ynom*H(1,1);
      gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
      PPLUS=(I-gain*H)*PMINUS;
      PPLUS=.5*(PPLUS+PPLUS');
      p11(i)=PPLUS(1,1);
      p22(i)=PPLUS(2,2);
      p12(i)=PPLUS(1,2);
      posvar(i)=p11(i)+2*ynom*p12(i)+ynom*ynom*p22(i);
      PMINUS=PHI*PPLUS*PHI'+Q;
   end
end

%Now do the plots.  Do bias, scale factor, and position separately.

k=0:1:2000;
plot(k,sqrt(p11))
title('Press ENTER to Continue')
pause

plot(k,sqrt(p22))
title('Press ENTER to Continue')
pause

plot(k,sqrt(posvar))
title('Press ENTER to End Problem')

%Comment: Note that the position error keeps increasing during
%the no-measurements period, even though the bias and scale factor
%errors are stable during this period.  The increase in position
%error is due to the multiplicative factor y on the scale factor
%error.  y continues to increase during the no-measurement period.

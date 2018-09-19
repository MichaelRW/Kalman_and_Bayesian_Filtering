clear
%File name is pb05_19.m and this is for Problem 5.19.

%First, specify the filter parameters.  Parameters are scalars.

sigma=1
beta=.25
dt=.02
phi=exp(-beta*dt)
q=sigma*sigma*(1-exp(-2*beta*dt))
sigmawk=sqrt(q)
h=1
r=1
pzerom=1

%Next, generate the 51-point z and true x sequences (4 of each).
%XTRUE and Z will be stacked 4 x 51 matrices.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0.  This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 26 of the code and replacing
%rand with randn everywhere.

rand('normal')
Z=zeros(4,51);
XTRUE=zeros(4,51);
for i=1:4
   XTRUE(i,1)=rand;
   Z(i,1)=XTRUE(i,1)+rand;
   for j=1:50
      XTRUE(i,j+1)=phi*XTRUE(i,j)+sigmawk*rand;
      Z(i,j+1)=XTRUE(i,j+1)+rand;
   end
end

%Next, do the filter estimates and error variances and stack the rows.

XHAT=zeros(4,51);
PBIG=zeros(4,51);
for i=1:4
   i
   xhatmin=0;
   pminus=pzerom;
   for j=1:51
      gain=pminus*h/(h*pminus*h+r);
      xhat=xhatmin+gain*(Z(i,j)-h*xhatmin);
      pplus=(1-gain*h)*pminus;
      XHAT(i,j)=xhat;
      PBIG(i,j)=pplus;
      xhatmin=phi*xhat;
      pminus=phi*pplus*phi+q;
   end
end

%To do the plots, we need both the positive and negative error sigmas
%and the Monte Carlo error sequences for the four estimates.

sigmapos=sqrt(PBIG(1,:));
sigmaneg=-sigmapos;
error1=XTRUE(1,:)-XHAT(1,:);
error2=XTRUE(2,:)-XHAT(2,:);
error3=XTRUE(3,:)-XHAT(3,:);
error4=XTRUE(4,:)-XHAT(4,:);
t=0:.02:1;
plot(t,sigmapos,'ow',t,sigmaneg,'ow',t,error1,t,error2,t,error3,t,error4)
title('Press ENTER to end Problem 5.19')

%To see more sample realizations, simply re-run pb05_19.m with
%a different seed.

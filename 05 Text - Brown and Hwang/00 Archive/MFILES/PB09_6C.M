clear
%File name is pb09_6c.m and this is for the R-bumped-up part
%of Problem 9.6.

%We will use Monte Carlo methods to estimate the mean square error.
%We will make many trial runs with the R-bumped-up filter and average
%the results.  The number of runs will be "runs" and it will be set
%at 1000, tentatively.  Set the number of steps s and runs:

s=51
runs=1000

%We will first precompute the gain sequence for the R-bumped-up filter.

PHIX=1
PHIMARK=exp(-.1)
QX=1
QMARK=1-exp(-.2)
H=1
R=2
PZEROM=1.618034
gain=zeros(1,s);
PMINUS=PZEROM;
for i=1:51
   K=PMINUS*H/(H*PMINUS*H+R);
   gain(i)=K;
   PPLUS=(1-K*H)*PMINUS;
   PMINUS=PHIX*PPLUS*PHIX+QX;
end

%Now make the Monte Carlo runs and sum the squared errors on a
%running basis to avoid a memory problem.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0.  This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 52 of the code and replacing
%rand with randn everywhere.

xtrue=zeros(1,s);
mark=zeros(1,s);
z=zeros(1,s);
xhatp=zeros(1,s);
xhatm=zeros(1,s);
error=zeros(1,s);
sqerror=zeros(1,s);
cumsqerr=zeros(1,s);
scale1=sqrt(1.618034);
scale2=sqrt(QMARK);
for i=1:runs
   rand('normal')
   i
   %Establish initial conditions at i=1 (i.e., t=0)
   xtrue(1)=scale1*rand;
   mark(1)=rand;
   z(1)=xtrue(1)+mark(1)+rand;
   xhatm(1)=0;
   xhatp(1)=(1-gain(1))*xhatm(1)+gain(1)*z(1);
   error(1)=xtrue(1)-xhatp(1);
   %Now do the remaining (s-1) steps recursively.
   for j=1:50
      xtrue(j+1)=PHIX*xtrue(j)+rand;
      mark(j+1)=PHIMARK*mark(j)+scale2*rand;
      z(j+1)=xtrue(j+1)+mark(j+1)+rand;
      xhatm(j+1)=PHIX*xhatp(j);
      xhatp(j+1)=(1-gain(j+1))*xhatm(j+1)+gain(j+1)*z(j+1);
      error(j+1)=xtrue(j+1)-xhatp(j+1);
   end

   %Now square the errors and accumulate squared errors in cumsqerr.
   
   sqerror=error.*error;
   cumsqerr=cumsqerr+sqerror;
end

%Finally, form the sequence of average squared errors.

avesqerr=(1/runs)*cumsqerr;

%Now bring p11opt and p11sk back into the workspace for plotting.
%Note that you must have run pb09_6a.m and pb09_6b.m (to generate
%the mat-files for p11opt and p11sk) to do this part.

load prob96a
load prob96b
t=0:1:(s-1);
plot(t,p11opt,'r',t,p11sk,'g',t,avesqerr,'w')
title('Press ENTER to end Problem 9.6')

%Note that the R-bumped-up filter performance is very comparable
%to the Schmidt-Kalman filter in this case.  This is because the
%Markov measurement "bias" is not a true bias and looks somewhat
%like white noise.  If it were truly white noise, then the
%R-bumped-up filter would, in fact, be optimum.

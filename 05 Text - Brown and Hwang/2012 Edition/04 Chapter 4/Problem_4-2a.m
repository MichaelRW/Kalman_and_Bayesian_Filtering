%This file is for Problem 4.2a
%Prob. to check triple integration model
%Specify paramters
W=1.0;
dt=10;   %single large step size
PHI=[1 dt (dt^2)/2;0 1 dt;0 0 1];
Q=[(1/20)*(dt^5) (1/8)*(dt^4) (1/6)*(dt^3);...
   (1/8)*(dt^4) (1/3)*(dt^3) (1/2)*(dt^2);...
   (1/6)*(dt^3) (1/2)*(dt^2) dt];
PMIN=zeros(3,3);
%No meas at t=0
PPLUS=PMIN;
PMIN10=PHI*PPLUS*PHI'+Q;
%Now repeat for 200 small steps
dt=.1;    %small incremental steps
PHI=[1 dt (dt^2)/2;0 1 dt;0 0 1];
Q=[(1/20)*(dt^5) (1/8)*(dt^4) (1/6)*(dt^3);...
   (1/8)*(dt^4) (1/3)*(dt^3) (1/2)*(dt^2);...
   (1/6)*(dt^3) (1/2)*(dt^2) dt];
PMIN=zeros(3,3);
pmin11=zeros(200,1);
pmin22=zeros(200,1);
pmin33=zeros(200,1);
%No meas at any of the smallsteps
for i=1:200
  PPLUS=PMIN;
  PMIN=PHI*PPLUS*PHI'+Q;
  pmin11(i)=PMIN(1,1);
  pmin22(i)=PMIN(2,2);
  pmin33(i)=PMIN(3,3);
end

  

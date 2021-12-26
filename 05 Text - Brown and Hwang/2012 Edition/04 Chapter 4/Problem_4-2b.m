%This file is for Problem 4.2b
%This is a continuation of Problem 4.2a
%Here we make a very crude approx. of Q
%Do only the 200 step part.  Specify parameters
W=1.0;
dt=.1;
PHI=[1 dt (dt^2)/2;0 1 dt;0 0 1];
Q=W*[0 0 0;0 0 0;0 0 dt];
PMIN=zeros(3,3);
pmin11=zeros(200,1);
pmin22=zeros(200,1);
pmin33=zeros(200,1);
for i=1:200
  PPLUS=PMIN;
  PMIN=PHI*PPLUS*PHI'+Q;
  pmin11(i)=PMIN(1,1);
  pmin22(i)=PMIN(2,2);
  pmin33(i)=PMIN(3,3);
end


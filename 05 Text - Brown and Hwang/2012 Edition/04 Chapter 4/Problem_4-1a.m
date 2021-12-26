%This file is for Problem 4.1a
%This is 2-state model
%State1 is bias, state2 is Wiener
PHI=[1 0;0 1];
Q=[0 0;0 1];
h=[1 1];
r=4;
PPLUS=zeros(2,2);
PMINUS=zeros(2,2);
var2=zeros(51,1);
PMINUS0=[100 0;0 0];
PMINUS=PMINUS0;
for i=1:51
  gain=PMINUS*h'*inv(h*PMINUS*h'+r);
  PPLUS=(eye(2)-gain*h)*PMINUS;
  %Now save variance of sum of x1 and x2
  var2(i)=PPLUS(1,1)+2*PPLUS(1,2)+PPLUS(2,2);
  %Project ahead
  PMINUS=PHI*PPLUS*PHI'+Q;
end

  

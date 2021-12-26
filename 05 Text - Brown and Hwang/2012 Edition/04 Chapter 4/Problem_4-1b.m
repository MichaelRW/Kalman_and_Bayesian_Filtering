%This file is for Problem 4.1b
%This is for 1-state model
phi=1.0;
q=1.0;
h=1.0;
r=4.0;
pminus0=100;
pplus=0.0;
var1=zeros(51,1);
pminus=pminus0;
for i=1:51
  gain=pminus*h*inv(h*pminus*h+r);
  pplus=(1.0-gain*h)*pminus;
  %Now save var1
  var1(i)=pplus;
  %Project ahead
  pminus=phi*pplus*phi+q;
end







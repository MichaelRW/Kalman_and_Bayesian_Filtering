%This file is for Problem 4.11
%Specify parameters
phi=1;
q=4;
r=1;
h=.5;    %for meas model 1
pminus0=100;
%Do meas model 1 first,then model 2
%Do update at t=0
gain0=pminus0*h/(h*pminus0*h+r);
pplus0=(1-gain0*h)*pminus0;
pminus=phi*pplus0*phi+q;
%Now do the remaining 200 steps
p1=zeros(200,1);
for i=1:200
  gain=pminus*h/(h*pminus*h+r);
  pplus=(1-gain*h)*pminus;
  %save pplus
  p1(i)=pplus;
  pminus=phi*pplus*phi+q;
end
%Next do meas 2 scenario
p2=zeros(200,1);
%Do update at t=0
gain20=pminus0*cos(1)/(cos(1)*pminus0*cos(1)+r);
pplus20=(1-gain20*cos(1))*pminus0;
pminus=phi*pplus20*phi+q;
%Now do remaining 200 steps for meas 2 model
for j=1:200
  h2=cos(1+(j/120));
  gain2=pminus*h2/(h2*pminus*h2+r);
  pplus2=(1-gain2*h2)*pminus;
  %save pplus2
  p2(j)=pplus2;
  pminus=phi*pplus2*phi+q;
end


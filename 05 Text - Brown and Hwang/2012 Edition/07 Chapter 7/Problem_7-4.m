clear
% This file is for Problem 7.4.
% This is a linearized covariance analysis problem and we are only
% interested in the diagonal terms of P.  The program will be similar
% to fildiag.m except H is time variable, so it will have to be
% computed inside the "for" loop.  Diagonal elements of P will be put
% into a column vector at each step and then stacked into a 4 x 201
% matrix called PDIAG.
% The key parameters have the same variable names as in fildiag.m.
% Now put key parameters in the workspace (all except H).

n=4
m=2
s=201
PHI=[1 0 0 0;0 1 0 0;0 0 exp(-.01) 0;0 0 0 exp(-.01)]
Q=[400 0 0 0;0 400 0 0;0 0 900*(1-exp(-.02)) 0;0 0 0 900*(1-exp(-.02))]
R=[225 0;0 225]
PZEROM=[2000 0 0 0;0 2000 0 0;0 0 900 0;0 0 0 900]

% Now set up the "for" loop for s steps and intialize certain variables.

I=eye(n);
PDIAG=zeros(n,s);
PMINUS=PZEROM;
d=10000;
for i=1:s
   i
   % Preliminary problem is to compute the time variable elements of H
   ynom=-10000+100*(i-1);
   rnom=sqrt(d^2+ynom^2);
   h11=d/rnom;
   h12=ynom/rnom;
   h21=-d/rnom;
   h22=ynom/rnom;
   H=[h11 h12 1 0;h21 h22 0 1];
   gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-gain*H)*PMINUS;
   PPLUS=.5*(PPLUS+PPLUS');
   for j=1:n
      PDIAG(j,i)=PPLUS(j,j);
   end
   PMINUS=PHI*PPLUS*PHI'+Q;
end

% Now make the plots.  Red is for x and green is for y.
t=0:1:200;
plot(t,PDIAG(1,:),'r',t,PDIAG(2,:),'g')
title('Press ENTER to end Problem 9.4')

% The y error covariance does not go to infinity when the nominal
% trajectory goes through the origin, because the assumed aircraft
% dynamics (i.e., random walk) gives the filter some memory, and all
% the weight is not placed on the current measurement.

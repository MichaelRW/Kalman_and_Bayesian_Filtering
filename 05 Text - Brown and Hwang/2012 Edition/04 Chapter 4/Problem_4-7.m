%This file is for Problem 4.7
%Specify the process parameters and use the vanLoan method
sigma=1.0;
omega0=.1;
dt=1.0;
W=1;
bsq=2*(sqrt(2))*omega0^3;
csq=(sqrt(2))*omega0^2;
bc=sqrt(bsq*csq);
scale=1.0+sqrt(2);
F=[0 1 0;0 0 1;-(omega0^3) -(scale)*(omega0^2) -(scale)*omega0];
G=[0;0;bc*sigma];
GWGT=[0 0 0;0 0 0;0 0 bsq*csq*W];
%Now do the vanLoan computation
A=dt*[-F GWGT;zeros(3,3) F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT';
Q=PHI*B(1:3,4:6);
C=(chol(Q))';
%Now specify initial process covariance
posvar=sigma^2;
velvar=posvar*(omega0^2)*(1.0/scale);
accvar=posvar*(omega0^4);
COVAR0=[posvar 0 0;0 velvar 0;0 0 accvar];
%Now compute a typical initial x vector
rand('normal')
rand('seed',1)
xzero=sigma*[rand;omega0*rand/(sqrt(scale));rand*(omega0^2)];
%Now generate the remaining samples of x recursively
X=zeros(3,101);
X(:,1)=xzero;
for i=1:100
  i
  X(:,i+1)=PHI*X(:,i)+C*[rand;rand;rand];
end
%Now plot the third row of X.  It should be well behaved.


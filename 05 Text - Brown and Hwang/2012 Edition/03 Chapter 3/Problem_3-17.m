%This file is for Problem 3.17
%Specify process parameters
sigma =1.0;
omega0=0.1;
zeta=1.0/sqrt(2);
num=(sigma^2)*(omega0^3)*2*sqrt(2);
dt=1.0;
%Now use van Loan method to get PHI and Q
%See Section 3.9
F=[0 1;-(omega0^2) -(2*zeta*omega0)];
G=[0;sqrt(num)];
W=1;
GWGT=G*G';
A=dt*[-F GWGT;zeros(2,2) F'];
B=expm(A);
PHI=B(3:4,3:4)';
Q=PHI*B(1:2,3:4);
%Now generate the wk sequence for t=1,t=2, ... t=100
%See Section 3.10
wk=zeros(2,100);
rand('normal')
rand('seed',1)
%We initialize x to be zero at t=0
x=zeros(2,101);
C=(chol(Q))';
for i=1:100
  wk(:,i)=C*[rand;rand];
end
%Now generate the desired x process
for j=1:100
  x(:,j+1)=PHI*x(:,j)+wk(:,j);
end



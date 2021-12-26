%This file is for Problem 7.6, Part (b)

clear all
 
randn('seed',12);           % Choose a randomization seed
xstate=[100; 0.3];          % State vector
Pcov=[45 0.2; 0.2 0.001];   % Error covariance matrix
R=[1 0; 0 1];               % Measurement noise matrix
 
N=2;                        % 2 state problem
alpha=1;
beta=2;
kappa=3-N;
lambda=alpha^2*(N+kappa)-N;
 
sqrtP=chol(Pcov)';          % Obtain Cholesky factor to generate random samples
 
xMC=zeros(2,5);
zMCpred=zeros(2,5);
sum=zeros(2,1);
L=5;
sigma=sqrt(lambda+N);
xMC(:,1)=xstate;                    % Computing the sigma points
xMC(:,2)=xstate+sqrtP(:,1)*sigma;
xMC(:,3)=xstate+sqrtP(:,2)*sigma;
xMC(:,4)=xstate-sqrtP(:,1)*sigma;
xMC(:,5)=xstate-sqrtP(:,2)*sigma;
for i=1:L
    zMCpred(1,i)=xMC(1,i)*cos(xMC(2,i));    % Generate predicted measurements
    zMCpred(2,i)=xMC(1,i)*sin(xMC(2,i));    %   corresponding to sigma points
    if i == 1
        sum=sum+lambda/(lambda+N).*zMCpred(:,i);
    else
        sum=sum+1/2/(lambda+N).*zMCpred(:,i);
    end;
end;
    
zpred=sum;
sumzzsq=zeros(2,2);
sumxzsq=zeros(2,2);
for i=1:L
    if i == 1
        sumxzsq=sumxzsq+(lambda/(lambda+N)+1-alpha^2+beta).*(xMC(:,i)-xstate)*(zMCpred(:,i)-zpred)';
        sumzzsq=sumzzsq+(lambda/(lambda+N)+1-alpha^2+beta).*(zMCpred(:,i)-zpred)*(zMCpred(:,i)-zpred)';
    else
        sumxzsq=sumxzsq+(1/2/(lambda+N)).*(xMC(:,i)-xstate)*(zMCpred(:,i)-zpred)';
        sumzzsq=sumzzsq+(1/2/(lambda+N)).*(zMCpred(:,i)-zpred)*(zMCpred(:,i)-zpred)';
    end;
end;
Pzzcov=sumzzsq+R;       % Eq. 7.5.13
Pxzcov=sumxzsq;         % Eq. 7.5.14
gain=Pxzcov/Pzzcov;     % Eq. 7.5.15 
disp(Pzzcov);
disp(gain);

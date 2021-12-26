%This file is for Problem 7.6, Part (a)

clear all
 
randn('seed',12);               % Choose a randomization seed
N=10000;                        % Total number of samples
xstate=[100; 0.3];              % State vector
Pcov=[45 0.2; 0.2 0.001];       % Error covariance matrix
R=[1 0; 0 1];                   % Measurement noise matrix
 
sqrtP=chol(Pcov)';              % Obtain Cholesky factor to generate random samples
 
xMC=zeros(2,N);                 % Initialize matrices
zMCpred=zeros(2,N);
sum=zeros(2,1);
sumsq=zeros(2,2);
for i=1:N
    xMC(:,i)=xstate+sqrtP*randn(2,1);       % Generate random state samples
    zMCpred(1,i)=xMC(1,i)*cos(xMC(2,i));    % Generate predicted measurements
    zMCpred(2,i)=xMC(1,i)*sin(xMC(2,i));    
    sum=sum+zMCpred(:,i);
    sumsq=sumsq+zMCpred(:,i)*zMCpred(:,i)';
end;
    
zpred=sum./N;
Pzzcov=sumsq./N-zpred*zpred'+R;     % Eq. 7.4.2
sumxzsq=zeros(2,2);
for i=1:N
    sumxzsq=sumxzsq+(xMC(:,i)-xstate)*(zMCpred(:,i)-zpred)';
end;
Pxzcov=sumxzsq./N;
gain=Pxzcov/Pzzcov;                 % Eq. 7.4.4
disp(Pzzcov);
disp(gain);

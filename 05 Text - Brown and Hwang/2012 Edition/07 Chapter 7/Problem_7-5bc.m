%This file is for Problem 7.5, Parts (b) and (c)

clear all
 
randn('seed',12);  
N=10000;                        % Total number of samples
xstate=[-16.5; -6.6; 2.9];      % State vector x
Pcov=[3.2 -1.2 0.2; -1.2 4.5 0.02; 0.2 0.02 1.2];   % Error covariance mtx
 
STM=[1 1 0;0 1 0;0 0 0.98];     % State transition matrix
Q=[1/3 1/2 0;1/2 1 0;0 0 0.01]; % Process noise matrix
sqrtQ=chol(Q)';
 
sqrtP=chol(Pcov)';
 
xMC=zeros(3,N);             % Initialize state vector
wMC=zeros(3,N);             % Initialize process noise vector
xMCproj=zeros(3,N);         % Initialize projected state vector
sum=0;
sumsq=0;
for i=1:N                                   % Monte Carlo loop
    xMC(:,i)=xstate+sqrtP*randn(3,1);       % Perturbing state vector  
    wMC(:,i)=sqrtQ*randn(3,1);              % Generating noise vector
    
    xMCproj(:,i)=STM*xMC(:,i)+wMC(:,i);     % Projected state vector sample
    sum=sum+xMCproj(:,i);                   % Summing sample
    sumsq=sumsq+xMCproj(:,i)*xMCproj(:,i)'; % Summing sample squared
end;
    
xproj=sum./N;                    % Computing mean of projected state samples
Pcovproj=sumsq./N-xproj*xproj';  % Computing covariance of projected state samples
 
disp(xproj);            % Display results
disp(Pcovproj);

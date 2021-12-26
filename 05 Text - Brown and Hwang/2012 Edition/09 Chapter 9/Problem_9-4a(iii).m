%This file is for Problem 9.4, Part (a), sub-part (iii)

%% Initialization
clear all;
 
randn('seed',2);        % Set the seed of randomization
dt=0.1;                 % Define sampling period
sigma=1;                % Sigma for the 1st order GM process
beta=1/3;               % Beta for the 1st order GM process
 
% State transition matrix and process noise covariance given by Ex. 3.12
phi=[1 (1-exp(-beta*dt))/beta; 0 exp(-beta*dt)];
q11=2*sigma^2/beta*(dt-2/beta*(1-exp(-beta*dt))+1/2/beta*(1-exp(-2*beta*dt)));
q12=2*sigma^2*(1/beta*(1-exp(-beta*dt))-1/2/beta*(1-exp(-2*beta*dt)));
q22=sigma^2*(1-exp(-2*beta*dt));
Q=[q11 q12; q12 q22];
A=chol(Q)';
 
%% Random process generation
N=100000;           % Total number of samples
x=zeros(2,N);
xx=zeros(1,N);
x(1,1)=sigma*randn(1,1);    % Initial states
x(2,1)=sigma*beta*randn(1,1);
bias=0.004;
for k=2:N
    x(:,k)=phi*x(:,k-1)+A*randn(2,1);
end;
xx=x(1,:);
 
%% Allan deviation
tau0=0.1;
for steps=0:15
    tau(steps+1)=tau0*2^steps;
    m=tau(steps+1)/tau0;
    sum=0;
    for i=1:N-2*m
        sum=sum+(xx(i+2*m)-2*xx(i+m)+xx(i))^2;
    end;
    sqrtav(steps+1)=sqrt(sum/2/(N-2*m)/tau(steps+1)^2);
end;
 
%% Plotting results
figure(1);
loglog(tau,sqrtav,'*-');
grid;
axis([0.1 1000 1e-3 1e1]);
xlabel('Averaging Time (s)');
ylabel('Allan deviation');

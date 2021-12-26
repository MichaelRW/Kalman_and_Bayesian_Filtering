%This file is for Problem 9.4, Part (a), sub-part (i)

%% Initialization
clear all;
randn('seed',2);  % Set the seed of randomization
 
dt=0.1;           % Define sampling period
sigma=1;          % Sigma for the 1st order GM process
sigmasq=sigma^2;
beta=1/3;         % Beta for the 1st order GM process
phi=exp(-beta*dt); % Evaluate the state transition term
qk=sigmasq*(1.0-exp(-2*beta*dt)); % Evaluate the process noise variance
C=sqrt(qk);
 
%% Random process generation
%Now form Monte Carlo xtrue sequence
%x1 and x2 are pos and vel phase variables
N=100000;          % Total number of samples
xtrue=zeros(2,N);  % Initial states
xtrue(1,1)=sigma*randn(1,1);
for i=2:N
    xtrue(1,i)=phi*xtrue(1,i-1)+C*randn(1,1);
end;
xx=xtrue(1,:);
 
%% Allan deviation
tau0=0.1;
for steps=0:14
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

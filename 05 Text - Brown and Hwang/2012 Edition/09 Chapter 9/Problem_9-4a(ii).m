%This file is for Problem 9.4, Part (a), sub-part (ii)

%% Initialization
clear all;
 
randn('seed',2);           % Set the seed of randomization
dt=0.1;                    % Define sampling period
sigma=1.0;                 % Sigma for the 2nd order GM process
sigmasq=sigma^2;
omega0=1/3;                % Natural freq. for the 2nd order GM process
bsq=2*(sqrt(2))*omega0^3;
b=sqrt(bsq);
 
%Compute state model PHI and Q
F=[0 1;-omega0^2 -(sqrt(2))*omega0]; % Forming F matrix for diff. eq.
GWGT=[0 0;0 bsq*sigmasq];            
A=dt*[-F GWGT;zeros(2,2) F'];        % Van Loan method for deriving
B=expm(A);                           % discrete-time solution parameters
PHI=B(3:4,3:4)';
Q=PHI*B(1:2,3:4);
 
%% Random process generation
%Now form Monte Carlo xtrue sequence
%x1 and x2 are pos and vel phase variables
C=chol(Q)';
N=100000;            % Total number of samples
xtrue=zeros(2,N);    % Initial states
xtrue(1,1)=sigma*randn(1,1);
xtrue(2,1)=sigma*omega0*randn(1,1);
for i=2:N
    xtrue(:,i)=PHI*xtrue(:,i-1)+C*randn(2,1);
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

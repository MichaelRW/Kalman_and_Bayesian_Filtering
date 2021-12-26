%This file is for Problem 9.3, Part (a)

clear all;
randn('seed',2);
dt=0.1;
phi=[1 dt; 0 1];
Sf=1E-21;   % Spectral amplitude of white noise driving input f
Sg=1.2E-22; % Spectral amplitude of white noise driving input g
Q=[Sf*dt+Sg/3*dt^3 Sg/2*dt^2; Sg/2*dt^2 Sg*dt];
C=chol(Q)';
N=100000;   % Total number of samples
x=zeros(2,N); % Initializing x vector
x(2,1)=1E-9;  % Initial state values
for k=2:N
    x(:,k)=phi*x(:,k-1)+C*randn(2,1);
end;
plot(x(1,:));

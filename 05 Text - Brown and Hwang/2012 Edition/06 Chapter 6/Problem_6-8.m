%This file is for Problem 6.8
%Discrete Wiener filter example
%For Example 6.5 and Problem 6.8
%Specify signal and noise parameters

sigmax=4.0;
sigmaxsq=sigmax^2;
sigmansq=2.0;
W=10;    %noise bandwidth
wx=1.0;  %signal bandwidth
psd=1/10;   %noise psd
dt=.05;   %sample spacing
n=50;     %number of samples
w=zeros(n,1);   %weighting vector
T=zeros(n,n);
f=zeros(n,1);
%Now form matrix parameters for solution
c=zeros(n,1);    %vector defining Toeplitz matrix
c(1)=sigmaxsq+sigmansq;
for i=1:(n-1)
  c(i+1)=sigmaxsq*(sin(2*pi*wx*i*dt))/(2*pi*wx*i*dt);
end
for i=1:(n-1)
  f(i)=sigmaxsq*(sin(2*pi*wx*(n-i)*dt))/(2*pi*wx*(n-i)*dt);
end
f(n)=sigmaxsq;
T=toeplitz(c);
%Now solve for w
w=(inv(T))*f;

%Finally, solve for mean square error
msqerror=sigmaxsq-2*w'*f+w'*T*w;


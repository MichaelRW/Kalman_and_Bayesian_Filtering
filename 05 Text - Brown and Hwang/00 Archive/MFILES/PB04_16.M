clear
%File name is pb04_16.m and this is for Problem 4.16.

%Analysis for this problem is similar to that for Problem 4.15.
%Autocorrelation function of s+n is exp(-tau) + exp(-4*tau).
%First form the Toeplitz matrix EZZ416.

ac416=zeros(32,1);
for i=1:32
   ac416(i)=exp(-.1*(i-1))+exp(-.1*4*(i-1));
end
EZZ416=toeplitz(ac416);

%Next, form the crosscorrelation vector ezs416.

ezs416=zeros(32,1);
for i=1:32
   ezs416(i)=exp(-.1*(32-i));
end

%Now solve for the weight factor k.

k=(inv(EZZ416))*ezs416;
k

%The result is consistent with the corresponding continuous Wiener
%problem where the weighting function is the sum of an impulse and
%an exponential.

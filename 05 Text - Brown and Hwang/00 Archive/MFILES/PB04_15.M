clear
%File name is pb04_15.m and this file is for Problem 4.15.

%We first form the z autocorrelation function vector from which
%the Toeplitz matrix EZZ is formed.

acsample=zeros(32,1);
for i=1:32
   acsample(i)=2*exp(-.1*(i-1));
end

%Now form EZZ.

EZZ=toeplitz(acsample);

%Next, form the zs crosscorrelation vector ezs.

ezs=zeros(32,1);
for i=1:32
   ezs(i)=exp(-.1*(32-i));
end

%Now solve for the weight factor vector k.

k=(inv(EZZ))*ezs;
k

%The resulting k is consistent with the corresponding continuous 
%Wiener solution where the weighting function is the impulse
%function .5*delta(t).

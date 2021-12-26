

%% Environment

close all force; clc; clear; restoredefaultpath;

set( 0, 'DefaultFigureWindowStyle', 'docked' );



%% ph

%This file is for Problem 2.25.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0.  This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 21 of the code and replacing
%rand with randn everywhere.

%Part (a):  First work out variance of wk and phi.

varwk=1-exp(-2*.05)
sigmawk=sqrt(varwk)
phi=exp(-.05)

%Set up a 1024-point sample realization of the process and call
%it sample (sample will be a 1 x 1024 row vector).

sample=zeros(1,1024);
rand( 0 );
sample(1)=rand;
for i=1:1023
   sample(i+1)=phi*sample(i)+sigmawk*rand;
end

%Now preview the sample realization for reasonableness.

plot(sample)


%Part(b):  First set up the true R(tau) for plotting.
%(Note that the first element of vector rtrue is for zero lag.)

tau=0:.05:3;
rtrue=zeros(1,61);
rtrue(1)=1;
for i=1:60
   rtrue(i+1)=exp(-.05*i);
end

%Next, we want to compute the experimental autocorrelation function.
%We do not wish to remove the sample mean, so we will write a short
%program to compute vxtau (called Vsubx(tau) in the text).

%s is the dimension of the row vector of time samples.
%m is the number of lags desired (not including zero lag).

s=1024
m=60
vxtau=zeros(1,m+1);
slag=zeros(1,s);
vxtau(1)=(sample*sample')/s;

%Note that the first element of vxtau is for zero lag.
%Now compute the remaining values of vxtau for s lags.

for i=1:m
   slag=[zeros(1,i) sample(1:(s-i))];
   vxtau(i+1)=(sample*slag')/(s-i);
end

%Now plot rtrue and vxtau together and compare.
        
plot(tau,vxtau,tau,rtrue)

%Part (c):  We will only do the ensemble of 8 vxtau's here and
%make a comparison of rtrue and the average of 8 vxtau's.  We
%already have one sample worked out.  We need 7 more.  They
%will be generated using similar code.  sample and vxtau will
%be re-used as dummy variables here.  The 8 vxtau's will be
%stacked in an 8 x 61 matrix called VXTAU.

VXTAU=zeros(8,61);
VXTAU(1,:)=vxtau;
for i=1:7
   i
   sample(1)=rand;
   for j=1:(s-1)
      sample(j+1)=phi*sample(j)+sigmawk*rand;
   end

   vxtau(1)=(sample*sample')/s;
   for j=1:m
      slag=[zeros(1,j) sample(1:(s-j))];
      vxtau(j+1)=(sample*slag')/(s-j);
   end

   VXTAU((i+1),:)=vxtau;
end

%Now plot all 8 vxtau's together and note the diversity.

plot(tau,VXTAU(1,:),tau,VXTAU(2,:),tau,VXTAU(3,:),tau,VXTAU(4,:),...
     tau,VXTAU(5,:),tau,VXTAU(6,:),tau,VXTAU(7,:),tau,VXTAU(8,:))
 

%Finally, form the average of the 8 vxtau's and compare the
%average with rtrue.  This is a vivid demonstration (not proof) 
%that the expectation of the experimental autocorrelation 
%function is equal to the true autocorrelation function.

vxtauav=mean(VXTAU);
plot(tau,vxtauav,tau,rtrue)



%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );





%% Environment

close all force; clc; clear; restoredefaultpath;

set( 0, 'DefaultFigureWindowStyle', 'docked' );



%% ph

%This file is for Problem 2.26.
                                                                    
%Parts (a) and (b):  We will first form 8 sample realizations of a
%Wiener process and plot them together for a reasonableness check.
%Then we will look at the average of the squares of 50 realizations
%and check to see if this is approximately linear with time.  The
%8 realizations can be conveniently stacked in an 8 x 11 matrix.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0.  This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 21 of the code and replacing
%rand with randn everywhere.

WIENER=zeros(8,11);

%Now generate the Wiener processes (note the initial value is zero).

randn( 0 );
for i=1:8
   for j=1:10
      WIENER(i,j+1)=WIENER(i,j)+rand;
   end
end

%Plot the first 8 samples and check for reasonableness.

t=0:1:10;
plot(t,WIENER(1,:),t,WIENER(2,:),t,WIENER(3,:),t,WIENER(4,:),...
     t,WIENER(5,:),t,WIENER(6,:),t,WIENER(7,:),t,WIENER(8,:))


%Part (c):  In this part we want to form the average of an 
%ensemble of squared Wiener-process realizations.  We will do
%the summing "on the run" and not save the identity of the
%individual runs.  The first 8 runs will be re-used, however.

%First square and sum the first 8 realizations.

SQUARED8=WIENER.*WIENER;
subsum8=sum(SQUARED8);

%Now add in the squared realizations for 42 more runs.

wiener=zeros(1,11);
parsum=subsum8;
for i=1:42
   for j=1:10
      wiener(j+1)=wiener(j)+rand;
   end
   wienersq=wiener.*wiener;
   parsum=parsum+wienersq;
end
mnsq50=(1/50)*parsum;
plot(t,mnsq50)




%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



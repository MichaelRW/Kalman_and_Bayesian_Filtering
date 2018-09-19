clear
%File name is pb03_26.m and this is for Problem 3.26.

%Part (a):  This part is analytical.  Pole is at .9512.
%Time constant is 1 sec, or 1/.05 = 20 steps.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0.  This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 27 of the code and replacing
%rand with randn everywhere.

%Part (b):  Similar to Prob. 2.33.  First work out phi and var(wk).

varwk=1-exp(-2*.05)
sigmawk= sqrt(varwk)
phi=exp(-.05)

%In this problem we need 4 sample realizations, each of length 501.
%Stack the samples as a 4 x 501 matrix.  (If this exceeds the matrix
%size limitation of your version of MATLAB, do each realization
%separately with different variable names for each.  Then modify
%the remaining code appropriately.)

SAMPLE=zeros(4,501);
rand('normal')

%The first element of each sample realization is set at zero.
%Now generate the remaining elements of SAMPLE.

for i=1:4
   for k=1:500
      SAMPLE(i,k+1)=phi*SAMPLE(i,k)+sigmawk*rand;
   end
end

%Plot the 4 SAMPLE row vectors together and check for reasonableness.

t=0:.05:25;
plot(t,SAMPLE(1,:),t,SAMPLE(2,:),t,SAMPLE(3,:),t,SAMPLE(4,:))
title('Press ENTER to Continue')
pause

%Part (c):  Start over and generate a new set of samples for the
%Wiener process.

varwk2=.095163
sigmawk2=sqrt(varwk2)
phi2=1

%Now generate the (4 x 501) matrix SAMPLE2 for part (c);

SAMPLE2=zeros(4,501);
for i=1:4
   for k=1:500
      SAMPLE2(i,k+1)=phi2*SAMPLE2(i,k)+sigmawk2*rand;
   end
end

%Now plot the results.

plot(t,SAMPLE2(1,:),t,SAMPLE2(2,:),t,SAMPLE2(3,:),t,SAMPLE2(4,:))
title('Press ENTER to end Problem 3.26')

%Note that the variance grows linearly with time, and thus the
%standard deviation increases according to the square root of 
%time.  Thus, the instability observed in the plots is not quite
%as pronounced as one might expect at first glance.

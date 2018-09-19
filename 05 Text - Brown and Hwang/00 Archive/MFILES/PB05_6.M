clear
%File name is pb05_6.m and this is for Problem 5.6.

%The purpose of this exercise is to provide a simple introduction to
%doing error covariance analysis using the Kalman filter recusive
%equations.  Some extra insight into detailed operation of a Kalman
%filter is gained by cycling through the recursive procedure one step
%at a time and observing the result with each step.  The result will
%be displayed as a plot that increments with each recursive step.

%First enter the key filter parameters into the workspace.

phi=exp(-.02)
q=1-exp(-2*.02)
r=1
h=1
pzerom=1

%Note that this example is scalar, so the recursive operations are
%scalar.  Set up the recursive loop to run 26 steps with a pause in
%the loop so the updated P can be observed as a plot with each step.
%(Also note that i=1 corresponds to t=0, i=2 to t=.02, and so forth.)

pvector=-ones(1,26);  %This initial pvector will not show on the plot

pminus=pzerom;
t=0:.02:.5;
for i=1:26
   gain=pminus*h/(h*pminus*h+r);
   pplus=(1-gain*h)*pminus;
   pvector(i)=pplus;
   axis([0 .5 0 .5]);
   plot(t,pvector,'o')
   axis([0 .5 0 .5]);
   if i<=25
      title('Press ENTER to Do Another Step')
   end
   if i==26
      title('Last step. Hit ENTER (and close plot window) for numerical data')
   end
   pause
   pminus=phi*pplus*phi+q;
end
save prob5_6.dat pvector /ascii
format compact
pvector
disp('    ')
disp('Reminder: Reset axis scaling to automatic before proceeding')
disp('to another problem.')

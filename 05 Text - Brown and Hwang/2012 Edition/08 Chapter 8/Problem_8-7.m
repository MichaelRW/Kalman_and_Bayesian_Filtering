clear

%This file is for Problem 8.7
%Part (a):  This part is analytic and the Q1 matrix is written out
%explicitly as a function of dt in the following code.  The PSD of
%the input white noise was set at unity for convenience.
 
%Parts (b) & (c):  We will now calculate both the Q1 and Q2 matrices
%for dt=5, dt=50, and dt=500.  First enter the key parameters.
 
g=9.8
re=6.37e6
F=[0 1 0;0 0 -g;0 1/re 0];
GWGT=[0 0 0;0 1 0;0 0 0];
 
%Now work out Q1 and Q2 for dt=5.
dt=5
Q1_5s=[(1/3)*dt*dt*dt       (1/2)*dt*dt     (1/3)*(1/re)*dt*dt*dt
       (1/2)*dt*dt               dt            (1/2)*(1/re)*dt*dt
    (1/3)*(1/re)*dt*dt*dt (1/2)*(1/re)*dt*dt (1/(3*re*re))*dt*dt*dt];
A=[-dt*F dt*GWGT;zeros(3) dt*F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT';
Q2_5s=PHI*B(1:3,4:6);
 
%Now work out Q1 and Q2 for dt=50.
dt=50
Q1_50s=[(1/3)*dt*dt*dt       (1/2)*dt*dt     (1/3)*(1/re)*dt*dt*dt
       (1/2)*dt*dt               dt            (1/2)*(1/re)*dt*dt
    (1/3)*(1/re)*dt*dt*dt (1/2)*(1/re)*dt*dt (1/(3*re*re))*dt*dt*dt];
A=[-dt*F dt*GWGT;zeros(3) dt*F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT';
Q2_50s=PHI*B(1:3,4:6);
 
%Now work out Q1 and Q2 for dt=500.
dt=500
Q1_500s=[(1/3)*dt*dt*dt       (1/2)*dt*dt     (1/3)*(1/re)*dt*dt*dt
       (1/2)*dt*dt               dt            (1/2)*(1/re)*dt*dt
    (1/3)*(1/re)*dt*dt*dt (1/2)*(1/re)*dt*dt (1/(3*re*re))*dt*dt*dt];
A=[-dt*F dt*GWGT;zeros(3) dt*F'];
B=expm(A);
PHIT=B(4:6,4:6);
PHI=PHIT';
Q2_500s=PHI*B(1:3,4:6);
 
%To avoid unnecessary clutter in our comparison, we will only
%compare corresponding diagonal terms in Q1 and Q2.  Also,
%because of large differences in units among the 11, 22, and
%33 terms, the comparisons for each will be made separately.
%We will make a diary of the results for perusal later.
 
s='1st row is Q1; 2nd row is Q2; columns are for dt=5, 50, & 500. '
format long
diary prob10_5
 
%First, compare the 11 terms (position variances).
 
comment=[s,'Comparison is for Q(1,1) terms.  Press ENTER to continue.']
Q11TERMS=[Q1_5s(1,1) Q1_50s(1,1) Q1_500s(1,1)
          Q2_5s(1,1) Q2_50s(1,1) Q2_500s(1,1)]
pause
 
%Next, compare the 22 terms.
 
comment=[s,'Comparison is for Q(2,2) terms.  Press ENTER to continue.']
Q22TERMS=[Q1_5s(2,2) Q1_50s(2,2) Q1_500s(2,2)
          Q2_5s(2,2) Q2_50s(2,2) Q2_500s(2,2)]
pause
 
%Finally, compare the 33 terms.
 
comment=[s,'Comparison is for Q(3,3) terms.  Press ENTER to complete prob.']
Q33TERMS=[Q1_5s(3,3) Q1_50s(3,3) Q1_500s(3,3)
          Q2_5s(3,3) Q2_50s(3,3) Q2_500s(3,3)]
pause
diary off
 
%Conclusion;  The approximation obtained by using the first-order
%phi in the Q calculation is excellent for dt=5 sec.  This is to be
%expected, because dt is only about 1/1000th of a Schuler period.
%The results are still quite good when dt is extended to 50 sec
%(1/100th of a Schuler period).  However, when we extend dt to
%500 sec (1/10th of a Schuler period), the accuracy degrades to
%somewhere in the 5 to 10 percent range which is poor by most any
%standards.

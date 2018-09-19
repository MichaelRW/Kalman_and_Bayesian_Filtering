clear
%File name is pb05_4.m and this is for Problem 5.4.

%Part (a):  This part is analytic.  The answers for F, GWG, PHI,
%and Q are given in Part (b).
 
%Part (b);  First, we need to compute PHI and Q for the model.
%Note that the white noise scale factor is accounted for in the
%output equation.

omega0=.012
beta=omega0/sqrt(2)
F=[0 1;-omega0^2 -omega0*sqrt(2)];
GWGT=[0 0;0 1];
A=[-10*F 10*GWGT;zeros(2,2) 10*F'];
B=expm(A);
PHIT=B(3:4,3:4);
PHI=PHIT'
Q=PHI*B(1:2,3:4)

%Now generate 4 sample realizations and stack the output y sequences.
%(If the stack of sample realizations exceeds the matrix size
%limitation of your MATLAB version, do only 2 realizations and
%change the subsequent code accordingly.)  Equations for the
%variances of x1 and x2 are found analytically.

%The rand statement used here is from MATLAB Version 3.5, and it
%results in a warning statement when using Version 4.0.  This does
%not interfere with the solution.  If you are using Version 4.0,
%the warning statement can be eliminated by deleting the
%rand('normal') statement in line 40 of the code and replacing
%rand with randn everywhere.

varx1=1/(2*sqrt(2)*omega0^3)
varx2=varx1*omega0^2
Y=zeros(4,401);
X=zeros(2,401);
CT=chol(Q);
C=CT'
rand('normal')
for i=1:4
   X(1,1)=(sqrt(varx1))*rand;
   X(2,1)=(sqrt(varx2))*rand;
   for j=1:400
      X(:,j+1)=PHI*X(:,j)+C*[rand;rand];
   end
   Y(i,:)=(sqrt(.002585))*X(1,:);
end
t=0:10:4000;
plot(t,Y(1,:),t,Y(2,:),t,Y(3,:),t,Y(4,:))
title('Press ENTER to Continue')
pause

%Part (c):  The code for finding the autocorrelation functions will
%be similar to that used in pb02_33.m.
                   
%s is the dimension of the vector of time samples.
%m is the number of lags (not including zero lag).

s=401
m=40
AUTOY=zeros(4,m+1);
vxtau=zeros(1,m+1);
ylag=zeros(1,s);
for i=1:4
   i
   vxtau(1)=(Y(i,:)*Y(i,:)')/s;

   %Note that the first element of vxtau is for zero lag.
   %Now compute the remaining elements of vxtau for m lags.

   for j=1:m
      ylag=[zeros(1,j) Y(i,1:(s-j))];
      vxtau(j+1)=(Y(i,:)*ylag')/(s-j);
   end
   AUTOY(i,:)=vxtau;
end

tau=0:10:400;
plot(tau,AUTOY(1,:),tau,AUTOY(2,:),tau,AUTOY(3,:),tau,AUTOY(4,:))
title('Press ENTER to Continue')
pause

%Part (d):  Finally, plot the average of the experimental auto-
%correlation functions along with the true a.c. function.

rtrue=zeros(1,41);
for i=1:41
   arg=beta*10*(i-1);
   rtrue(i)=23*23*(exp(-arg))*(cos(arg)+sin(arg));
end
aveautoy=mean(AUTOY);
plot(tau,aveautoy,tau,rtrue,'o')
title('Press ENTER to end Problem 5.4')

%If you want to see more sample realizations, simply repeat
%pb05_4.m using a different seed.

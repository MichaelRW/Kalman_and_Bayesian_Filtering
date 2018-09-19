clear
%File name is pb06_7.m and this is for Problem 6.7.

%We will use the methodology given in the first part of Section 6.7
%and analyze the three filter scenarios shown in Fig. 6.10.  This
%may not be the most efficient way to program the problem, but it
%it is easier to understand this way.

%First enter key parameters for the optimal filter (truth model).

s=96
omega=2*pi*60
dt=1/(60*32)
omegadt=pi/16
PZEROM=[.5 0 0 0 0 0;0 .5 0 0 0 0;0 0 .005 0 0 0;0 0 0 .005 0 0;...
        0 0 0 0 .00125 0;0 0 0 0 0 .00125];
PHI=eye(6);
Q=[1/3200 0 0 0 0 0;0 1/3200 0 0 0 0;0 0 1/320000 0 0 0;...
   0 0 0 1/320000 0 0;0 0 0 0 1/1280000 0;0 0 0 0 0 1/1280000];
R=[.0009]
I=eye(6);

%H varies with each step, so it will be specified inside the "for" loop.
%We will only save the diagonal elements of the optimal P sequence.
%The diagonal elements will be stacked side-by-side as vectors in OPTP.

OPTP=zeros(6,s);
PMINUS=PZEROM;
for i=1:s
   k=i-1
   H=[cos(k*omegadt) -sin(k*omegadt) cos(3*k*omegadt) -sin(3*k*omegadt) ...
      cos(5*k*omegadt) -sin(5*k*omegadt)];
   GAIN=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-GAIN*H)*PMINUS;
   %Symmetrize PPLUS
   PPLUS=(PPLUS+PPLUS')/2;
   for j=1:6
      OPTP(j,i)=PPLUS(j,j);
   end
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Now run the suboptimal filter (Filter No. 2 in Fig. 6.10) with the
%wrong Q.  The suboptimal gains will be saved in SUBGAIN.

SUBGAIN=zeros(6,s);
PMINUS=PZEROM;
Q(3,3)=4*Q(3,3);
Q(4,4)=4*Q(4,4);
for i=1:s
   k=i-1
   H=[cos(k*omegadt) -sin(k*omegadt) cos(3*k*omegadt) -sin(3*k*omegadt) ...
      cos(5*k*omegadt) -sin(5*k*omegadt)];
   GAIN=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-GAIN*H)*PMINUS;
   %Symmetrize PPLUS
   PPLUS=(PPLUS+PPLUS')/2;
   SUBGAIN(:,i)=GAIN;
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Finally, we re-cycle the suboptimal gains through the truth
%model (Filter No. 1 in Fig. 6.10).  We will stack the diagonal
%elements of the resulting P in a matrix called SUBP.

Q(3,3)=.25*Q(3,3);
Q(4,4)=.25*Q(4,4);
SUBP=zeros(6,s);
PMINUS=PZEROM;
for i=1:s
   k=i-1
   H=[cos(k*omegadt) -sin(k*omegadt) cos(3*k*omegadt) -sin(3*k*omegadt) ...
      cos(5*k*omegadt) -sin(5*k*omegadt)];
   GAIN=SUBGAIN(:,i);
   PPLUS=(I-GAIN*H)*PMINUS*(I-GAIN*H)'+GAIN*R*GAIN';
   for j=1:6
      SUBP(j,i)=PPLUS(j,j);
   end
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Now plot OPTP and SUBP together to make a comparison.  Ignore the
%transient during the first cycle and plot data for the last 2 cycles.

t=32:1:(s-1);
plot(t,OPTP(1,33:s),'r',t,SUBP(1,33:s),'g')
title('Press ENTER to Continue')
pause
plot(t,OPTP(2,33:s),'r',t,SUBP(2,33:s),'g')
title('Press ENTER to Continue')
pause
plot(t,OPTP(3,33:s),'r',t,SUBP(3,33:s),'g')
title('Press ENTER to Continue')
pause
plot(t,OPTP(4,33:s),'r',t,SUBP(4,33:s),'g')
title('Press ENTER to Continue')
pause
plot(t,OPTP(5,33:s),'r',t,SUBP(5,33:s),'g')
title('Press ENTER to Continue')
pause
plot(t,OPTP(6,33:s),'r',t,SUBP(6,33:s),'g')
title('Press ENTER to end Problem 6.7')

%The 4 to 1 variation in Q(3,3) and Q(4,4) does not affect the 
%estimation error in the fundamental and fifth harmonic appreciably.
%The effect in the third harmonic is about 20 percent on the basis
%of variance, or about 10 percent on an rms basis.  This is a
%relatively low sensitivity situation.

clear
%File name is pb06_5.m and this is for Problem 6.5.

%We only need the gain sequence for this problem, so we will not
%need to generate a sample realization of the process as was done
%in Section 6.5.

%First enter key parameters.  s is the number of steps.

s=160
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

%H varies with each step, so it will be specified inside the "for"
%loop.  We will save the gain sequence in a stacked matrix GAINSEQ.

GAINSEQ=zeros(6,s);
PMINUS=PZEROM;
for i=1:s
   k=i-1
   H=[cos(k*omegadt) -sin(k*omegadt) cos(3*k*omegadt) -sin(3*k*omegadt) ...
      cos(5*k*omegadt) -sin(5*k*omegadt)];
   GAIN=(PMINUS*H')*inv(H*PMINUS*H'+R);
   PPLUS=(I-GAIN*H)*PMINUS;
   %Symmetrize PPLUS
   PPLUS=(PPLUS+PPLUS')/2;
   GAINSEQ(:,i)=GAIN;
   PMINUS=PHI*PPLUS*PHI'+Q;
end

%Plot gain sequence for fundamental.

t=0:1:(s-1);
plot(t,GAINSEQ(1,:),'r',t,GAINSEQ(2,:),'g')
title('Press ENTER to Continue')
pause

%Plot gain sequence for third harmonic.

plot(t,GAINSEQ(3,:),'r',t,GAINSEQ(4,:),'g')
title('Press ENTER to Continue')
pause

%Plot gain sequence for fifth harmonic.

plot(t,GAINSEQ(5,:),'r',t,GAINSEQ(6,:),'g')
title('Press ENTER to end Problem 6.5')

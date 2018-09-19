clear
%File name is pb10_2.m and this is for Problem 10.2.

%First enter key parameters into workspace and work out PHI and Q.
%Begin with the no-accel. case & then do the -4 m/sec/sec next.

n=9
m=2
s=201
dt=1
g=9.8
re=6.37e6
omegax=1.57e-5
omegay=7.27e-5
wacc=.0036
wgyro=2.35e-9
R=[225 0;0 225]
ay=4.0
PZEROM=zeros(n);
PZEROM(1,1)=100;
PZEROM(2,2)=1e-6;
PZEROM(3,3)=1e-6;
PZEROM(4,4)=100;
PZEROM(5,5)=1e-6;
PZEROM(6,6)=1e-6;
PZEROM(7,7)=100;
PZEROM(8,8)=1e-6;
PZEROM(9,9)=(1.0*pi/180)^2

%Now work on PHI and Q matrices.  Start with zero acceleration.

F=zeros(n);
F(1,2)=1;
F(2,3)=-g;
F(3,2)=1/re;
F(3,9)=omegax;
F(4,5)=1;
F(5,6)=-g;
F(6,5)=1/re;
F(6,9)=omegay;
F(7,8)=1;

GWGT=zeros(n);
GWGT(2,2)=wacc;
GWGT(3,3)=wgyro;
GWGT(5,5)=wacc;
GWGT(6,6)=wgyro;
GWGT(8,8)=wacc;
GWGT(9,9)=wgyro;

A=[-dt*F dt*GWGT;zeros(n) dt*F'];
B=expm(A);
PHIT=B(10:18,10:18);
PHI=PHIT';
Q=PHI*B(1:9,10:18);

%Now order constant negative accel. steps k=95,96,...104 (in k.F.
%notation beginning with k=0), as i=96,97,...105 in the "for" loop.
%Call phi and q in mid-region PHIMID and QMID. Now calculate these.

F(2,9)=-ay
AMID=[-dt*F dt*GWGT;zeros(n) dt*F'];
BMID=expm(AMID);
PHIMIDT=BMID(10:18,10:18);
PHIMID=PHIMIDT';
QMID=PHIMID*BMID(1:9,10:18);

%Begin the main i loop.  PDIAG will contain the diagonal
%elements of PPLUS in a 9 x 201 stacked matrix.  (If this exceeds
%the matrix size limitation of the version of MATLAB that you are
%using, try running the dynamic scenario from 50 sec to 150 sec.
%This will yield nearly the same results as for the whole 200 sec.)

PMINUS=PZEROM;
PDIAG=zeros(9,201);
H=zeros(2,9);
d=10000
I=eye(n);

for i=1:201
   i
   %Compute variable H matrix.  Linearize about a trajectory
   %that is symmetric about the midpoint (i.e., the origin).
   if i>0 & i<=95
      ynom=-10000+100*(i-1);
      rnom=sqrt(d^2+ynom^2);
      H(1,1)=d/rnom;
      H(1,4)=ynom/rnom;
      H(2,1)=-d/rnom;
      H(2,4)=ynom/rnom;
      gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
      PPLUS=(I-gain*H)*PMINUS;
      PPLUS=.5*(PPLUS+PPLUS');
      for j=1:9
         PDIAG(j,i)=PPLUS(j,j);
      end
      PMINUS=PHI*PPLUS*PHI'+Q;
   end

   if i>95 & i<=105
      ynom=-10000+80*(i-1);
      rnom=sqrt(d^2+ynom^2);
      H(1,1)=d/rnom;
      H(1,4)=ynom/rnom;
      H(2,1)=-d/rnom;
      H(2,4)=ynom/rnom;
      gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
      PPLUS=(I-gain*H)*PMINUS;
      PPLUS=.5*(PPLUS+PPLUS');
      for j=1:9
         PDIAG(j,i)=PPLUS(j,j);
      end
      PMINUS=PHIMID*PPLUS*PHIMID'+QMID;
   end

   if i>105 & i<=201
      ynom=-10000+60*(i-1);
      rnom=sqrt(d^2+ynom^2);
      H(1,1)=d/rnom;
      H(1,4)=ynom/rnom;
      H(2,1)=-d/rnom;
      H(2,4)=ynom/rnom;
      gain=(PMINUS*H')*inv(H*PMINUS*H'+R);
      PPLUS=(I-gain*H)*PMINUS;
      PPLUS=.5*(PPLUS+PPLUS');
      for j=1:9
         PDIAG(j,i)=PPLUS(j,j);
      end
      PMINUS=PHI*PPLUS*PHI'+Q;
   end
end

%Now plot azimuth estimation error variance

k=0:1:200;
plot(k,PDIAG(9,:))

%Comment:  Note the dramatic decrease in azimuth error during and
%immediately following the acceleration pulse.  The drop is 
%nearly 4 to 1 on a variance basis (2 to 1 rms), and this is for
%an acceleration pulse that is relatively modest both in amplitude
%and duration.  Perhaps even more important is the fast response,
%relative to what could be obtained with gyrocompassing.

%This is a general covariance analysis program for a Schmidt-
%Kalman filter.  The final result is the error covariance for x
%in the form of stacked matrix called COVARX.

%File name is sch_kal.m

%Declare the global variables.

global n r m s PHIX PHIY QX QY H J R PX0M PXY0M PY0M COVARX

%n is order of x; r is the order of y; m is the order of z;
%s is the number of steps (including k=0).

I=zeros(n,n);
for j=1:n
   I(j,j)=1
end

COVARX=zeros(n,n*s);

%post M stands for minus and post P stands for plus.
%Initialize P

PXM=PX0M;
PXYM=PXY0M;
PYXM=PXYM';
PYM=PY0M;
for i=1:s
   i
   ALPHA=H*PXM*H'+H*PXYM*J'+J*PYXM*H'+J*PYM*J'+R;
   KX=(PXM*H'+PXYM*J')*(inv(ALPHA));
   PXP=(I-KX*H)*PXM-KX*J*PYXM;
   PXYP=(I-KX*H)*PXYM-KX*J*PYM;
   PYXP=PXYP';
   PYP=PYM;
   COVARX(1:n,(n*(i-1)+1):n*i)=PXP;
   PXM=PHIX*PXP*PHIX'+QX;
   PXYM=PHIX*PXYP*PHIY';
   PYXM=PXYM';
   PYM=PHIY*PYP*PHIY'+QY;
end

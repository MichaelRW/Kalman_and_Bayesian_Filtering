clear
% File name is PB11_8A.M and this is for Problem 11.8 part (a).

% This m-file uses a function file, DRNCOS1.M that provides the direction
% cosine matrix for four satellites specified in Problem 11.8.
% First, specify the filter parameters.

  c=3e8;
  dt=1;
  TotalSteps=300;
  state_dim=5;
  PHI=eye(state_dim);
  PHI(4,5)=dt;
  Sp=0.01;
  Sf=0.4E-18*c^2;
  Sg=1.58E-18*c^2;
  Q=[Sp*dt  0       0         0               0;
     0    Sp*dt     0         0               0;
     0      0     Sp*dt       0               0;
     0      0       0  Sf*dt+(Sg*dt^3)/3  (Sg*dt^2)/2;
     0      0       0      (Sg*dt^2)/2      Sg*dt];

  R=0.01*eye(4);


  PMINUS=[10000 0 0 0 0;
          0 10000 0 0 0;
          0 0 10000 0 0;
          0 0 0 90000 0;
          0 0 0  0  900];

  Hfill=[1 0;1 0;1 0;1 0];
  drncos0=drncos1(0);      % function call to DRNCOS1.M.
  time=0:1:TotalSteps;
  Save_P=zeros(3,TotalSteps+1);

  disp ('Running...');

% Main recursive loop.
  for k=1:TotalSteps+1
    time(k)=k-1;
    H=[drncos1(k-1)-drncos0  Hfill];
    GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
    PPLUS=(eye(state_dim)-GAIN*H)*PMINUS;
    PPLUS=(PPLUS+PPLUS')/2;
    PMINUS=PHI*PPLUS*PHI'+Q;
    PMINUS=(PMINUS+PMINUS')/2;
    Save_P(1,k)=PPLUS(1,1);
    Save_P(2,k)=PPLUS(2,2);
    Save_P(3,k)=PPLUS(3,3);
  end;

%Scaling below is for MATLAB 4.  Scaling was avoided for MATLAB S3.5
%because of unresolved problems when involving semilog plots.
  semilogy(time,Save_P(1,:),time,Save_P(2,:),time,Save_P(3,:));
  axis([0 TotalSteps 1e-2 1e4]);

%Return graphics to auto scaling mode.
  clear axis;

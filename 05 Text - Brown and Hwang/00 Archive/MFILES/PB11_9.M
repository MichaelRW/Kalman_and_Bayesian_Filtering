clear
% File name is PB11_9.M and this is for Problem 11.9 part (a).

% This m-file uses a function file, DRNCOS2.M that provides the direction
% cosine matrix for four satellites specified in Problem 11.9.
% First, specify the filter parameters.

  c=3e8;
  dt=1;
  TotalSteps=1200;
  state_dim=4;
  PHI=eye(state_dim);
  Q=[1e10  0    0    0;
     0    1e10  0    0;
     0     0   1e10  0;
     0     0    0   1e10];

  R=eye(4);


  PMINUS=[1e10 0 0 0;
          0 1e10 0 0;
          0 0 1e10 0;
          0 0 0 1e10];

  time=200:1:TotalSteps;
  Save=zeros(1,TotalSteps-199);
  Hfill=[1;1;1;1];

  disp ('Running...');

%Starting index of 200 avoids limitation of MATLAB S3.5
  for k=200:TotalSteps
    H=[drncos2(k-1)  Hfill];
    GAIN=PMINUS*H'*inv(H*PMINUS*H'+R);
    PPLUS=(eye(state_dim)-GAIN*H)*PMINUS;
    PPLUS=(PPLUS+PPLUS')/2;
    Save(k-199)=sqrt(PPLUS(1,1)+PPLUS(2,2)+PPLUS(3,3)+PPLUS(4,4));
    PMINUS=PHI*PPLUS*PHI'+Q;
    PMINUS=(PMINUS+PMINUS')/2;
  end;

%Plot the GDOP over t=200 to 1200
  semilogy(time,Save);
  title('Press ENTER to end Problem 11.9');
  pause;

%Output to the workspace, the GDOP value at t=1200
  disp('GDOP=');
  disp(sqrt(PPLUS(1,1)+PPLUS(2,2)+PPLUS(3,3)+PPLUS(4,4)));
